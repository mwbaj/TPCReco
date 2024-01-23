#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <iomanip>
#include <unistd.h>
#include "gtest/gtest.h"
#include "TFile.h"

#include "TPCReco/EventSourceFactory.h"
#include "TPCReco/ConfigManager.h"
#include "TPCReco/colorText.h"

#include "TPCReco/grawToEventTPC.h"


class grawToEventTPCTest : public ::testing::Test {
public:
  static std::shared_ptr<PEventTPC> myEventPtr;

  static TFile* rootfile;

  static boost::property_tree::ptree myConfig;

  static void SetUpTestSuite() {
  
    // load the graw file
    std::string testJSON = std::string(std::getenv("HOME"))+"/.tpcreco/config/test.json";
    int argc = 3;
    char *argv[] = {(char*)"ConfigManager_tst", 
                  (char*)"--meta.configJson",const_cast<char *>(testJSON.data())};
                
    ConfigManager cm;
    myConfig = cm.getConfig(argc, argv);
    int status = chdir("../../resources");
    std::shared_ptr<EventSourceBase> myEventSource = EventSourceFactory::makeEventSourceObject(myConfig);
  
    myEventPtr = myEventSource->getCurrentPEvent(); 
    myEventSource->loadFileEntry(0);

    // convert the graw file to root file
    convertGRAWFile(myConfig);

    rootfile = new TFile("EventTPC_2022-04-12T08_03_44.531_0000.root", "READ");
  }

  static void TearDownTestSuite() {
    std::remove("EventTPC_2022-04-12T08_03_44.531_0000.root");
  }


};

std::shared_ptr<PEventTPC> grawToEventTPCTest::myEventPtr(0);
TFile* grawToEventTPCTest::rootfile(0);
boost::property_tree::ptree grawToEventTPCTest::myConfig;


TEST(ROOTFileNameTest, createROOTFileName)
{
  std::string grawFileName = "../testData/CoBo0_AsAd0_2022-04-12T08_03_44.531_0000.graw,"
                             "../testData/CoBo0_AsAd1_2022-04-12T08_03_44.533_0000.graw,"
                             "../testData/CoBo0_AsAd2_2022-04-12T08_03_44.536_0000.graw,"
                             "../testData/CoBo0_AsAd3_2022-04-12T08_03_44.540_0000.graw";
  std::string rootFileName = createROOTFileName(grawFileName);
  EXPECT_EQ(rootFileName, "EventTPC_2022-04-12T08_03_44.531_0000.root");
}

TEST_F(grawToEventTPCTest, fileExists) 
{
    ASSERT_EQ(rootfile->IsZombie(), false);
    EXPECT_EQ(rootfile->IsOpen(), true);
    rootfile->Close();
    EXPECT_EQ(rootfile->IsOpen(), false);
}

TEST_F(grawToEventTPCTest, compareGrawToRoot)
{
  auto grawChargeMap = myEventPtr->GetChargeMap();
  TTree* tree = (TTree*)rootfile->Get(myConfig.get<std::string>("input.treeName").c_str());
  
  std::shared_ptr<EventSourceBase> rootEventSource = EventSourceFactory::makeEventSourceObject(myConfig);
  std::shared_ptr<PEventTPC> rootEventPtr = rootEventSource->getCurrentPEvent();
  rootEventSource->loadFileEntry(0);
  auto rootChargeMap = rootEventPtr->GetChargeMap();
  EXPECT_EQ(grawChargeMap, rootChargeMap);
}

