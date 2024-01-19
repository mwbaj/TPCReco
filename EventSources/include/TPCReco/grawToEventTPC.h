#ifndef _grawToEventTPC_H_
#define _grawToEventTPC_H_

#include <iostream>
#include <cstdlib>
#include <string>
#include <memory>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/program_options.hpp>


#include "TPCReco/EventTPC.h"
#include "TPCReco/PEventTPC.h"
#include "TPCReco/EventSourceGRAW.h"
#include "TPCReco/EventSourceMultiGRAW.h"
#include "TPCReco/EventSourceFactory.h"
#include "TPCReco/ConfigManager.h"

#include <TFile.h>
#include <TTree.h>



std::string createROOTFileName(const  std::string & grawFileName);

int convertGRAWFile(boost::property_tree::ptree & aConfig);

#endif