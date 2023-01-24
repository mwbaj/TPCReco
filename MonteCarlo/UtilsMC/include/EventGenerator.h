#ifndef TPCSOFT_EVENTGENERATOR_H
#define TPCSOFT_EVENTGENERATOR_H

#include "PrimaryParticle.h"
#include "ReactionLibrary.h"
#include "ZProvider.h"
#include "XYProvider.h"
#include "GeneratorSetup.h"
#include "pugixml.h"
#include "boost/property_tree/ptree.hpp"

class EventGenerator {
public:
    EventGenerator()=delete;
    EventGenerator(boost::property_tree::ptree configNode);
    EventGenerator(const boost::filesystem::path &configFilePath);
private:
    void Init();
    ReactionLibrary lib{};
    std::unique_ptr<ZProvider> zProv;
    std::unique_ptr<XYProvider> xyProv;
    std::unique_ptr<EProvider> eProv;
    GeneratorSetup setup;
    ROOT::Math::XYZPoint beamPosition;
    double beamAngleZ{};
};


#endif //TPCSOFT_EVENTGENERATOR_H
