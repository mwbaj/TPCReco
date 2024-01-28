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
    // edit for reproducing fitted events from GUI - if the changes are accepted, possibly as a new module, the values should also be taken from the config file (example ones put here for now)
    TVector3 vertex = {-49.8369, 0.868339, -37.4044}; 
    TVector3 prod1End = {-25.522, -59.766, -54.9762};
    TVector3 prod2End = {-52.3815, 8.26315, -34.3788};
    double prod1length = 67.6498;
    double prod2length = 8.38527;
    TVector3 prod1tangent = {1, 1.83356, -1.18942};
    TVector3 prod2tangent = {1, 1.20164, 1.90221};
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
