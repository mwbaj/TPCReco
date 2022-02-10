#include "TrackSegment2D.h"
#include "colorText.h"
#include <iostream>

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void TrackSegment2D::setBiasTangent(const TVector3 & aBias, const TVector3 & aTangent){

  myBias = aBias;
  myTangent = aTangent.Unit();

  double lambda = 10;//FIXME what value should be here?
  myStart = myBias;
  myEnd = myBias+lambda*myTangent;
  initialize();
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void TrackSegment2D::setStartEnd(const TVector3 & aStart, const TVector3 & aEnd){

  myStart = aStart;
  myEnd = aEnd;

  myTangent = (myEnd - myStart).Unit();
  myBias = (myStart + myEnd)*0.5;

  initialize();
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void TrackSegment2D::initialize(){
  
  double lambda = -myBias.X()/myTangent.X();
  myBiasAtT0 = myBias + lambda*myTangent;
  
  lambda = -myBias.Y()/myTangent.Y();
  myBiasAtStrip0 = myBias + lambda*myTangent;

  ///Set tangent direction along time arrow with unit time component,
  ///so vector components can be compared between projections.
  myTangentWithT1 = myTangent;
  if(myTangentWithT1.X()<0) myTangentWithT1 *= -1;
  
  if(std::abs(myTangentWithT1.X())>1E-5){
    myTangentWithT1 *= 1.0/myTangentWithT1.X();
  }
  else{
    myTangentWithT1.SetZ(0.0);
    myTangentWithT1 = myTangentWithT1.Unit();
  }
  
  myLenght = (myEnd - myStart).Mag();
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
const TGraph & TrackSegment2D::getChargeProfile(const Hit2DCollection & aRecHits, double radiusCut){

  double x = 0.0, y = 0.0;
  double distance = 0.0;
  double lambda = 0.0;
  TVector3 aPoint;
  myChargeProfile.Set(0);
  
  for(const auto aHit:aRecHits){
    x = aHit.getPosTime();
    y = aHit.getPosStrip();
    aPoint.SetXYZ(x, y, 0.0);
    std::tie(lambda, distance) = getPointLambdaAndDistance(aPoint);
    if(lambda>0 && lambda<getLength() && distance<radiusCut){
      myChargeProfile.SetPoint(myChargeProfile.GetN(), lambda,  aHit.getCharge());
    }
  }
  return myChargeProfile;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
double TrackSegment2D::getIntegratedCharge(double lambdaCut, const Hit2DCollection & aRecHits) const{

  double x = 0.0, y = 0.0;
  double totalCharge = 0.0;
  double radiusCut = 4.0;//FIXME put into configuration. Value 4.0 abtained by  looking at the plots.
  double distance = 0.0;
  double lambda = 0.0;
  TVector3 aPoint;
  
  for(const auto aHit:aRecHits){
    x = aHit.getPosTime();
    y = aHit.getPosStrip();
    aPoint.SetXYZ(x, y, 0.0);    
    std::tie(lambda, distance) = getPointLambdaAndDistance(aPoint);
    if(lambda>0 && lambda<getLength() && distance>0 && distance<radiusCut) totalCharge += aHit.getCharge();
  }
  return totalCharge;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::tuple<double,double> TrackSegment2D::getPointLambdaAndDistance(const TVector3 & aPoint) const{

  const TVector3 & tangent = getTangent();
  const TVector3 & start = getStart();

  TVector3 delta = aPoint - start;
  double lambda = delta*tangent/tangent.Mag();
  TVector3 transverseComponent = delta - lambda*tangent;
  return std::make_tuple(lambda, transverseComponent.Mag());
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
double TrackSegment2D::getRecHitChi2(const Hit2DCollection & aRecHits) const {

  double dummyChi2 = 1E9;

  if(getTangent().Mag()<1E-3){
    std::cout<<__FUNCTION__<<KRED<< " TrackSegment2D has null tangent "<<RST
	     <<" for direction: "<<getStripDir()<<std::endl;
    std::cout<<KGRN<<"Start point: "<<RST;
    getStart().Print();
    std::cout<<KGRN<<"End point: "<<RST;
    getEnd().Print();
    return dummyChi2;
  }

  TVector3 aPoint;
  double chi2 = 0.0;
  double biasDistance = 0.0;
  double chargeSum = 0.0;
  double distance = 0.0;
  double lambda = 0.0;
  double minLambda = 0.0;
  double maxLambda = 0.0;
  double deltaLambda = 0.0;
  int pointCount = 0;
  
  double x = 0.0, y = 0.0;
  double charge = 0.0;

  for(const auto aHit:aRecHits){
    x = aHit.getPosTime();
    y = aHit.getPosStrip();
    charge = aHit.getCharge();
    aPoint.SetXYZ(x, y, 0.0);
    std::tie(lambda,distance) = getPointLambdaAndDistance(aPoint);
    if(lambda<0 || lambda>getLength()) continue;
    ++pointCount;
    chi2 += std::pow(distance, 2)*charge;
    chargeSum +=charge;
    biasDistance += (aPoint - getBias()).Mag();
    if(lambda<minLambda) minLambda = lambda;
    if(lambda>maxLambda) maxLambda = lambda;
  }
  if(!pointCount) return dummyChi2;

  chi2 /= chargeSum;
  biasDistance /= chargeSum;
  deltaLambda = std::abs(maxLambda - minLambda);
  deltaLambda = 0;
  //biasDistance = 0;
  return chi2 + biasDistance + deltaLambda;//TEST
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::ostream & operator << (std::ostream &out, const TrackSegment2D &aSegment){

  const TVector3 & start = aSegment.getStart();
  const TVector3 & end = aSegment.getEnd();
  
  out<<"direction: "<<aSegment.getStripDir()
     <<" ("<<start.X()<<", "<<start.Y()<<")"
     <<" -> "
     <<"("<<end.X()<<", "<<end.Y()<<") "<<std::endl
     <<"\t N Hough accumulator hits: "
     <<"["<<aSegment.getNAccumulatorHits()<<"]";
  return out;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
