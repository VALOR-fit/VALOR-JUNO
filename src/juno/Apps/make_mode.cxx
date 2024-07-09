#include "valor/Conventions/ReacMode.h"
#include "valor/GENIEUtilities/UtlEventRec.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <map>
#include <iostream>
#include <vector>

// Generates a fixed set of colors for the histograms
std::vector<int> generateDistinctColors() {
    return {kRed, kBlue, kGreen, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal};  // Predefined colors from ROOT
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root> <branch_name>\n";
        return 1;
    }

    const char* inputFileName = argv[1];
    const char* branchName = argv[2];

    TFile* file = TFile::Open(inputFileName);
    if (!file) {
        std::cerr << "Failed to open the file: " << inputFileName << "\n";
        return 1;
    }

    TTree* tree = nullptr;
    file->GetObject("valor_tree", tree);
    if (!tree) {
        std::cerr << "Failed to find the tree 'valor_tree' in " << inputFileName << "\n";
        file->Close();
        delete file;
        return 1;
    }

    if (tree->GetBranch(branchName) == nullptr) {
        std::cerr << "Branch '" << branchName << "' not found in the tree.\n";
        file->Close();
        delete file;
        return 1;
    }

    valor::conv::ReacMode_t mode;
    double value;
    tree->SetBranchAddress("reac_mode", &mode);
    tree->SetBranchAddress(branchName, &value);

    std::vector<int> colors = generateDistinctColors();
    std::map<valor::conv::ReacMode_t, TH1F*> histograms;
    
    for (int i = 0; i < valor::conv::kNReacModes; ++i) {
        valor::conv::ReacMode_t currentMode = static_cast<valor::conv::ReacMode_t>(i);
        std::string title = valor::util::AsString(currentMode);
    // Correcting histogram title setup to ensure it's meaningful and included in the legend
        histograms[currentMode] = new TH1F(title.c_str(), (title + ";Energy;Counts").c_str(), 5000, 0.0, 16.0);
        if (title == "nueCCoth") {
            histograms[currentMode]->SetFillColor(colors[0]);  // First color for nueCCoth
        } else {
            histograms[currentMode]->SetFillColor(colors[(i % (colors.size() - 1)) + 1]);  // Start from second color
        }
            histograms[currentMode]->GetXaxis()->SetRangeUser(0.0, 20.0);
    }

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (histograms.count(mode) > 0) {
            histograms[mode]->Fill(value);
        }
    }

THStack* stack = new THStack("stack", "Q^{2} Distributions for Different Reaction Modes");

// Set titles before drawing


TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
for (auto& hist : histograms) {
    if (hist.second->GetEntries() > 0) {
        stack->Add(hist.second);
        legend->AddEntry(hist.second, hist.second->GetTitle(), "f");
    }
}

TCanvas* canvas = new TCanvas("c1", "Canvas for Stacked Histogram", 800, 600);
canvas->SetLogy();
stack->Draw("bar"); // Draw the stack after setting all properties

legend->Draw();
canvas->Modified(); // Ensure the canvas is updated
canvas->Update(); // Additional update to ensure changes are visible

canvas->SaveAs("Q2_distributions.png");

file->Close();
delete file;
delete canvas;
delete legend;

return 0;
}
