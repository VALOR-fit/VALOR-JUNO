#ifndef JUNO_SINGLE_SAMPLE_FIT_MODEL_H
#define JUNO_SINGLE_SAMPLE_FIT_MODEL_H

#include "TXMLEngine.h"

#include "valor/FitElements/GenericSingleSampleFitModel.h"
#include "JunoDetectorResponse.h"

namespace valor {

class JunoSingleSampleFitModel : public GenericSingleSampleFitModel {
public:
  JunoSingleSampleFitModel();
  virtual ~JunoSingleSampleFitModel();

protected:
  bool ReadExperimentSpecificConfig(
      TXMLEngine* xml,
      XMLNodePointer_t root_node) override;

  bool ApplyExperimentSpecificResponse(
      Option_t* opt) override;

  bool UseDummyRecoForPreResponse() const override;

private:
  bool fUseJunoDetectorResponse;
  juno::JunoDetectorResponse fJunoDetectorResponse;
  
};

} // namespace valor

#endif