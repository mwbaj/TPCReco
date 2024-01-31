#include "TPCReco/IonRangeCalculator.h"
#include "TPCReco/ReactionTwoProng.h"
#include "Math/Boost.h"
#include "Math/LorentzRotation.h"
#include "TVector3.h"

using namespace std::string_literals;

ReactionTwoProng::ReactionTwoProng(std::unique_ptr<AngleProvider> theta, std::unique_ptr<AngleProvider> phi,
                                   pid_type targetIon, pid_type prod1Ion, pid_type prod2Ion)
        : thetaProv{std::move(theta)}, phiProv{std::move(phi)}, target{targetIon}, prod1Pid{prod1Ion},
          prod2Pid{prod2Ion} {
    CheckStoichiometry({targetIon},{prod1Ion,prod2Ion});
    targetMass = ionProp->GetAtomMass(target);
    prod1Mass = ionProp->GetAtomMass(prod1Ion);
    prod2Mass = ionProp->GetAtomMass(prod2Ion);
}



PrimaryParticles
ReactionTwoProng::GeneratePrimaries(double gammaMom, const ROOT::Math::Rotation3D &beamToDetRotation) {
    GetKinematics(gammaMom, targetMass);
    auto Qvalue = totalEnergy - prod1Mass - prod2Mass;
    //return empty vector if we do not have enough energy
    if (Qvalue < 0) {
        auto msg = "Beam energy is too low to create "s;
        msg += enumDict::GetPidName(prod1Pid) + "+" + enumDict::GetPidName(prod2Pid);
        msg += " pair! Creating empty event";
        std::cout << msg << std::endl;
        return {};
    }
    auto p_CM = 0.5 * sqrt((pow(totalEnergy, 2) - pow(prod2Mass - prod1Mass, 2)) *
                           (pow(totalEnergy, 2) - pow(prod1Mass + prod2Mass, 2))) / totalEnergy;
    auto thetaCM = thetaProv->GetAngle();
    auto phiCM = phiProv->GetAngle();
    // edit for reproducing fitted events from GUI - if changes are accepted, possibly as a new module, the values should also be taken from the config file (example ones put here for now)
    TVector3 vertex = {-150.013, 1.41374, -7.32393}; 
    TVector3 prod1End = {-108.383, 2.34073, -54.2624};
    TVector3 prod2End = {-155.324, 2.25, -0.956189};
    double prod1length = 62.7468;
    double prod2length = 8.33388;
    TVector3 prod1tangent = {1, 2.41593, 0.0222633};
    TVector3 prod2tangent = {1, 0.701183, 2.98542};
    IonRangeCalculator ionCalculator;
    double prod1_T_LAB = ionCalculator.getIonEnergyMeV(prod1Pid,prod1length);
    double prod2_T_LAB = ionCalculator.getIonEnergyMeV(prod2Pid,prod2length);

// leftover from previous version
    ROOT::Math::Polar3DVector p3(p_CM, thetaCM, phiCM);
    ROOT::Math::PxPyPzEVector p4FirstCM(p3.X(), p3.Y(), p3.Z(),
                                               sqrt(p_CM * p_CM + prod1Mass * prod1Mass));
    ROOT::Math::PxPyPzEVector p4SecondCM(-p3.X(), -p3.Y(), -p3.Z(),
                                                sqrt(p_CM * p_CM + prod2Mass * prod2Mass));

    ROOT::Math::Boost bst(-betaCM);
    ROOT::Math::LorentzRotation lRot(beamToDetRotation);
//
    PrimaryParticles result;

    double prod1_p_LAB=sqrt(prod1_T_LAB*(prod1_T_LAB+2*prod1Mass));
    double prod2_p_LAB=sqrt(prod2_T_LAB*(prod2_T_LAB+2*prod2Mass));

    TVector3 firstDir = (prod1End-vertex).Unit();
    TVector3 secondDir = (prod2End-vertex).Unit();

    TLorentzVector p4First(prod1_p_LAB*firstDir, prod1Mass+prod1_T_LAB);
    TLorentzVector p4Second(prod2_p_LAB*secondDir, prod2Mass+prod2_T_LAB);

    result.emplace_back(prod1Pid, p4First);
    result.emplace_back(prod2Pid, p4Second);
    return result;

}
