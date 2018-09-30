#include "cell.h"

const double CCell::vc=-80;
//const double CCell::stim=80;
//const double CCell::stimduration=2;

// ---------------constant parameters ------------------------------
const double CCell::temp=308.0;// temperature (K)
const double CCell::xxr=8.314;//
const double CCell::xf=96.485;// Faraday's constant
const double CCell::frt=xf/(xxr*temp);

CCell::CCell(void) : y(new double[N]),
m(y[0]), h(y[1]), j(y[2]), xr(y[3]),
xs1(y[4]), xs2(y[5]), xtos(y[6]), ytos(y[7]),
v(y[8]), cp(y[9]), cs(y[10]), ci(y[11]),
cnsr(y[12]), cjsr(y[13]), xir(y[14]), c1(y[15]),
c2(y[16]), xi1ca(y[17]), xi1ba(y[18]), xi2ca(y[19]),
xi2ba(y[20]), nai(y[21]), xtof(y[22]), ytof(y[23]),
tropi(y[24]), trops(y[25]) {

  // initial conditions
  m = 0.001145222753;// sodium m-gate
  h = 0.9898351676;// sodium h-gate
  j = 0.9930817518;// sodium j-gate

  xr = 0.008709989976;// ikr gate variable 
  xs1 = 0.08433669901;// iks gate variable
  xs2 = 0.1412866149;// iks gate varaible 

  xtos = 0.003757746357;// ito slow activation
  ytos = 0.1553336368;// ito slow inactivation

  v = -86.79545769; // voltage

  cp = 1.682601371;// averaged dyadic space con.
  cs = 0.3205609256;// averaged submembrane conc.
  ci = 0.3863687451;// myoplasm conc.

  cnsr = 107.0388739;// NSR load
  cjsr = 95.76256179;// average JSR load

  xir = 0.006462569526;// SR current flux

  // Markov gate variables 

  c1 = 1.925580885e-05;// C1
  c2 = 0.9535940241;// C2
  xi1ca = 0.007052299702;// I1_Ca
  xi1ba = 3.629261123e-05;// I1_Ba
  xi2ca = 0.02316349806;// I2_Ca
  xi2ba = 0.01613268649;// I2_Ba

  nai = 14.01807252;// internal Na conc.

  xtof = 0.003737842131;// ito fast activation
  ytof = 0.9823715315;// ito slow inactivation

  tropi = 29.64807803;// time dependent buffers in myplasm (troponin)
  trops = 26.37726416;// time dependent buffers in submembrane (troponin)

  dt = 0.1;
  vold = v;
  jparam = 1;

  nao = 136.0;//mM   external Na
  ki = 140.0;// mM  internal K
  ko = 5.40;//mM    external K
  cao = 1.8;// mM    external Ca

  ek = (1.0 / frt)*log(ko / ki);// K reversal potential


  gca = 182;// ica conductance
  gtos = 0.04;// ito slow conductance 
  gtof = 0.11;// ito fast conductance 
  gnaca = 0.84;// exchanger strength 
  gkr = 0.0125;// Ikr conductance 
  gks = 0.1386;
  gkix = 0.3;// Ik1 conductance
  gnak = 1.5;
  vup = 0.4;//0.3;// uptake strength
  taud = 4.0;// diffusional delay (ms)
  gna = 12.0;// sodium conductance (mS/micro F) 
  taur = 30.0;// spark lifetime (ms)
  taua = 100.0;// NSR-JSR diffusional delay (ms)
  av = 11.3;
  cstar = 90.0;
  gleak = 0.00002069;
  knsr = 50;
  cup = 0.5;// uptake threshold
}
CCell::~CCell() {
  delete[] y;
}
CCell& CCell::operator=(const CCell& cell)
{
  if (&cell != this) {
    for (int i = 0; i < N; i++) {
      y[i] = cell.y[i];
    }
    jparam = cell.jparam;
    vold = cell.vold;
    dt = cell.dt;

    nao = cell.nao;
    ki = cell.ki;
    ko = cell.ko;
    cao = cell.cao;

    gca = cell.gca;
    gtos = cell.gtos;
    gtof = cell.gtof;
    gnaca = cell.gnaca;
    gkr = cell.gkr;
    gks = cell.gks;
    gkix = cell.gkix;
    gnak = cell.gnak;
    vup = cell.vup;
    taud = cell.taud;
    gna = cell.gna;
    taur = cell.taur;
    taua = cell.taua;
    av = cell.av;
    cstar = cell.cstar;
    gleak = cell.gleak;
    knsr = cell.knsr;
    cup = cell.cup;
  }
  return(*this);
}

double CCell::pace(double st) {
  // -------------time step adjustment ------------------------
  double dv = (vold - v) / dt;
  vold = v;
  double itotal;
  if (fabs(dv) > 25.0) {// then finer time step when dv/dt large
    dtx = dt / 10;
    for (int iii = 0; iii < 10; iii++) {
      itotal = pacex(st);
    }
  }
  else {
    dtx = dt;
    itotal = pacex(st);
  }
  return itotal;
}


double CCell::pacex(double st) {
  ik1 = comp_ik1();
  double ito = comp_ito();//itos and itof
  inak = comp_inak();
  double csm = cs / 1000.0;// convert micro M to mM
  double jnaca = comp_jnaca(csm);
  //----------- Equations for Ca cycling -------------------------
  double jd = (cs - ci) / taud;//diffusion from submembrane to myoplasm
  // Troponin kinetics
  const double xkon = 0.0327;
  const double xkoff = 0.0196;
  const double btrop = 70.0;
  double xbi = xkon*ci*(btrop - tropi) - xkoff*tropi;
  double xbs = xkon*cs*(btrop - trops) - xkoff*trops;

  double jup = comp_jup();
  double jleak = comp_jleak();

  double po = comp_icalpo();
  double rxa = comp_rxa(csm);
  double jca = gca*po*rxa;// Ca current in micro M/ms

  double dcs = comp_inst_buffer(cs)*(50.0*(xir - jd - jca + jnaca) - xbs);
  double dci = comp_inst_buffer(ci)*(jd - jup + jleak - xbi);
  double dcnsr = -xir + jup - jleak;// SR load dynamics 
  double dcjsr = (cnsr - cjsr) / taua;// NSR-JSR relaxation dynamics 
  double Qr = comp_Q();
  double dir = comp_dir(po, Qr, rxa, dcnsr);
  double dcp = comp_dcp(po, Qr, rxa);

  ina = comp_ina();
  ikr = comp_ikr();
  iks = comp_iks();

  cp += dcp*dtx;
  cs += dcs*dtx;
  ci += dci*dtx;
  cnsr += dcnsr*dtx;
  xir += dir*dtx;
  cjsr += dcjsr*dtx;

  tropi += xbi*dtx;
  trops += xbs*dtx;

  //-------convert ion flow to current---------
  const double wca = 8.0;//conversion factor between micro molar/ms to micro amps/ micro farads
  inaca = wca*jnaca;
  ica = 2.0*wca*jca;
  //--------sodium dynamics -------------------------
  const double xrr = (1.0 / wca) / 1000.0;// note: sodium is in m molar so need to divide by 1000
  nai += (-xrr*(ina + 3.0*inak + 3.0*inaca))*dtx;
  // -------- dV/dt ------------------------------------
  double itotal = (-(ina + ik1 + ikr + iks + ito + inaca + ica + inak) + st);
  v += itotal*dtx;

  return itotal;
}
//-----------   sodium current following Hund-Rudy -------------------
double CCell::comp_ina(void) {
  double ena = (1.0 / frt)*log(nao / nai);
  double am = 3.2;
  if (fabs(v + 47.13) > 0.001)am = 0.32*(v + 47.13) / (1.0 - exp(-0.1*(v + 47.13)));
  double bm = 0.08*exp(-v / 11.0);
  double ah = 0;
  double aj = 0;
  double bh, bj;
  if (v < -40.0) {
    ah = 0.135*exp((80.0 + v) / (-6.8));
    bh = 3.56*exp(0.079*v) + 310000.0*exp(0.35*v);
    aj = ((-127140.0*exp(0.2444*v) - 0.00003474*exp(-0.04391*v))*(v + 37.78)) / (1.0 + exp(0.311*(v + 79.23)));
    bj = (0.1212*exp(-0.01052*v)) / (1.0 + exp(-0.1378*(v + 40.14)));
  }
  else {
    bh = 1.0 / (0.13*(1.0 + exp((v + 10.66) / (-11.1))));
    bj = (0.3*exp(-0.0000002535*v)) / (1.0 + exp(-0.1*(v + 32.0)));
  }

  double tauh = 1.0 / (ah + bh);
  double tauj = 1.0 / (aj + bj)*jparam;
  double taum = 1.0 / (am + bm);
  double xina = gna*h*j*m*m*m*(v - ena);

  h = ah / (ah + bh) - ((ah / (ah + bh)) - h)*exp(-dtx / tauh);
  j = aj / (aj + bj) - ((aj / (aj + bj)) - j)*exp(-dtx / tauj);
  m = am / (am + bm) - ((am / (am + bm)) - m)*exp(-dtx / taum);
  return xina;
}
//-------------- Ikr following Shannon------------------ 
double CCell::comp_ikr(void) {
  const double gss = sqrt(ko / 5.4);
  double xkrv1 = 0.00138 / 0.123;
  if (fabs(v + 7.0) > 0.001)
    xkrv1 = 0.00138*(v + 7.0) / (1. - exp(-0.123*(v + 7.0)));
  double xkrv2 = 0.00061 / 0.145;
  if (fabs(v + 10.0) > 0.001)
    xkrv2 = 0.00061*(v + 10.0) / (exp(0.145*(v + 10.0)) - 1.0);
  double taukr = 1.0 / (xkrv1 + xkrv2);
  double xkrinf = 1.0 / (1.0 + exp(-(v + 50.0) / 7.5));
  double rg = 1.0 / (1.0 + exp((v + 33.0) / 22.4));
  double xikr = gkr*gss*xr*rg*(v - ek);
  xr = xkrinf - (xkrinf - xr)*exp(-dtx / taukr);
  return xikr;
}
// ----- Iks modified from Shannon, with new Ca dependence------------
double CCell::comp_iks(void) {
  const double prnak = 0.018330;
  double eks = (1.0 / frt)*log((ko + prnak*nao) / (ki + prnak*nai));
  double xs1ss = 1.0 / (1.0 + exp(-(v - 1.50) / 16.70));
  double xs2ss = xs1ss;
  double tauxs1;
  if (fabs(v + 30.0) < 0.001 / 0.0687)
    tauxs1 = 1 / (0.0000719 / 0.148 + 0.000131 / 0.0687);
  else
    tauxs1 = 1.0 / (0.0000719*(v + 30.0) / (1.0 - exp(-0.148*(v + 30.0))) + 0.000131*(v + 30.0) / (exp(0.0687*(v + 30.0)) - 1.0));
  double tauxs2 = 4 * tauxs1;
  double gksx = (1 + 0.8 / (1 + pow((0.5 / ci), 3)));
  double xiks = gks*gksx*xs1*xs2*(v - eks);
  xs1 = xs1ss - (xs1ss - xs1)*exp(double(-dtx / tauxs1));
  xs2 = xs2ss - (xs2ss - xs2)*exp(double(-dtx / tauxs2));
  return xiks;
}
//------Ik1 following Luo-Rudy formulation (from Shannon model) ------
double CCell::comp_ik1(void) {
  const double gki = (sqrt(ko / 5.4));
  double aki = 1.02 / (1.0 + exp(0.2385*(v - ek - 59.215)));
  double bki = (0.49124*exp(0.08032*(v - ek + 5.476)) + exp(0.061750*(v - ek - 594.31))) / (1.0 + exp(-0.5143*(v - ek + 4.753)));
  double xkin = aki / (aki + bki);
  double xik1 = gkix*gki*xkin*(v - ek);
  return xik1;
}
//------- Ito slow following Shannon et. al. 2005 -----------
//------- Ito fast following Shannon et. al. 2005 -----------
double CCell::comp_ito(void) {
  double rt1 = -(v + 3.0) / 15.0;
  double rt2 = (v + 33.5) / 10.0;
  double rt3 = (v + 60.0) / 10.0;
  double xtos_inf = 1.0 / (1.0 + exp(rt1));
  double ytos_inf = 1.0 / (1.0 + exp(rt2));
  double rs_inf = 1.0 / (1.0 + exp(rt2));
  double txs = 9.0 / (1.0 + exp(-rt1)) + 0.5;
  double tys = 3000.0 / (1.0 + exp(rt3)) + 30.0;
  itos = gtos*xtos*(ytos + 0.5*rs_inf)*(v - ek);// ito slow
  xtos = xtos_inf - (xtos_inf - xtos)*exp(-dtx / txs);
  ytos = ytos_inf - (ytos_inf - ytos)*exp(-dtx / tys);

  double xtof_inf = xtos_inf;
  double ytof_inf = ytos_inf;
  double rt4 = -(v / 30.0)*(v / 30.0);
  double rt5 = (v + 33.5) / 10.0;
  double txf = 3.5*exp(rt4) + 1.5;
  double tyf = 20.0 / (1.0 + exp(rt5)) + 20.0;
  itof = gtof*xtof*ytof*(v - ek);// ito fast
  xtof = xtof_inf - (xtof_inf - xtof)*exp(-dtx / txf);
  ytof = ytof_inf - (ytof_inf - ytof)*exp(-dtx / tyf);
  return itos + itof;
}
// -------Inak (sodium-potassium exchanger) following Shannon --------------
double CCell::comp_inak(void) {
  const double xkmko=1.5; //these are Inak constants adjusted to fit
                          //the experimentally measured dynamic restitution curve
  const double xkmnai=12.0;
  const double sigma = (exp(nao/67.3)-1.0)/7.0;
  double fnak = 1.0/(1+0.1245*exp(-0.1*v*frt)+0.0365*sigma*exp(-v*frt));
  double xinak = gnak*fnak*(1./(1.+(xkmnai/nai)))*ko/(ko+xkmko);
  return xinak;
}
// --- Inaca (sodium-calcium exchange) following Shannon and Hund-Rudy------
//  Note: all concentrations are in mM
double CCell::comp_jnaca(double csm) {
  double zw3 = pow(nai, 3)*cao*exp(v*0.35*frt) - pow(nao, 3)*csm*exp(v*(0.35 - 1.)*frt);
  double zw4 = 1.0 + 0.2*exp(v*(0.35 - 1.0)*frt);
  const double xkdna = 0.3;// micro M
  double aloss = 1.0 / (1.0 + pow((xkdna / cs), 3));
  const double xmcao = 1.3;
  const double xmnao = 87.5;
  const double xmnai = 12.3;
  const double xmcai = 0.0036;
  double yz1 = xmcao*pow(nai, 3) + pow(xmnao, 3)*csm;
  double yz2 = pow(xmnai, 3)*cao*(1.0 + csm / xmcai);
  double yz3 = xmcai*pow(nao, 3)*(1.0 + pow((nai / xmnai), 3));
  double yz4 = pow(nai, 3)*cao + pow(nao, 3)*csm;
  double zw8 = yz1 + yz2 + yz3 + yz4;
  double xinacaq = gnaca*aloss*zw3 / (zw4*zw8);
  return xinacaq;
}
//  compute driving force
double CCell::comp_rxa(double csm) {
  const double pca = 0.00054;
  double za = v*2.0*frt;
  double factor1 = 4.0*pca*xf*xf / (xxr*temp);
  double factor = v*factor1;
  double rxa;
  if (fabs(za) < 0.001) {
    rxa = factor1*(csm*exp(za) - 0.341*(cao)) / (2.0*frt);
  }
  else {
    rxa = factor*(csm*exp(za) - 0.341*(cao)) / (exp(za) - 1.0);
  }
  return rxa;
}
// ------ Markovian Ca current --------------------------------
//  Markov model:All parameters have been fitted directly to 
//  experimental current traces using a multidimensional current fitting
//  routine.    
double CCell::comp_icalpo(void) {
  const double vth=0.0;
  const double s6=8.0;
  const double taupo=1.0;
  double poinf=1.0/(1.0+exp(-(v-vth)/s6));
  
  double alpha=poinf/taupo;
  double beta=(1.0-poinf)/taupo;

  const double r1=0.30;
  const double r2=3.0;
  const double cat=3.0;
  double fca=1.0/(1.0+pow(double(cat/cp),3));
  double s1=0.0182688*fca;
  const double s1t=0.00195;
  double k1=0.024168*fca;
  const double k2=1.03615e-4;
  const double k1t=0.00413;
  const double k2t=0.00224;
  double s2=s1*(r1/r2)*(k2/k1);
  const double s2t=s1t*(r1/r2)*(k2t/k1t);

  const double vx=-40;
  const double sx=3.0;
  double poi=1.0/(1.0+exp(-(v-vx)/sx));
  const double tau3=3.0;
   
  double k3=(1.0-poi)/tau3;
  double k3t=k3;
          
  const double vy=-40.0;
  const double sy=4.0;
  double Pr=1.0-1.0/(1.0+exp(-(v-vy)/sy));

  double recov=10.0+4954.0*exp(v/15.6);

  const double tca=78.0329;
  const double cpt=6.09365;
  double tau_ca=tca/(1.0+pow((cp/cpt),4))+0.1;

  double tauca=(recov-tau_ca)*Pr+tau_ca;
  double tauba=(recov-450.0)*Pr+450.0;

  const double vyr=-40.0;
  const double syr=11.32;
  double Ps=1.0/(1.0+exp(-(v-vyr)/syr));

  double k6=fca*Ps/tauca;
  double k5=(1.0-Ps)/tauca;
       
  double k6t=Ps/tauba;
  double k5t=(1.0-Ps)/tauba;

  double k4=k3*(alpha/beta)*(k1/k2)*(k5/k6);
  double k4t=k3t*(alpha/beta)*(k1t/k2t)*(k5t/k6t);

  double po=1.0-xi1ca-xi2ca-xi1ba-xi2ba-c1-c2;
  double dc2= beta*c1+k5*xi2ca+k5t*xi2ba-(k6+k6t+alpha)*c2;
  double dc1=alpha*c2+k2*xi1ca+k2t*xi1ba+r2*po-(beta+r1+k1t+k1)*c1;
  double dxi1ca=k1*c1+k4*xi2ca+s1*po-(k3+k2+s2)*xi1ca;
  double dxi2ca=k3*xi1ca+k6*c2-(k5+k4)*xi2ca;
  double dxi1ba=k1t*c1+k4t*xi2ba+s1t*po-(k3t+k2t+s2t)*xi1ba;
  double dxi2ba=k3t*xi1ba+k6t*c2-(k5t+k4t)*xi2ba;

  c1+=dc1*dtx;
  c2+=dc2*dtx;
  xi1ca+=dxi1ca*dtx;
  xi1ba+=dxi1ba*dtx;
  xi2ca+=dxi2ca*dtx;
  xi2ba+=dxi2ba*dtx;
  return po;
}
//----- SERCA2a uptake current ------------------------------------
double CCell::comp_jup(void) {
  return vup*ci*ci/(ci*ci + cup*cup);
}
// ---------leak from the SR--------------------------
double CCell::comp_jleak(void) {
  return gleak*(cnsr*cnsr / (cnsr*cnsr + knsr*knsr))*(cnsr*16.667 - ci);//vsr/vcell=0.06
}
// ---------- buffer dynamics in the myoplasm -----------------------
//buffering to calmodulin and SR are instantaneous, while buffering to
//Troponin C is time dependent.These are important to have reasonable
//Ca transient.Note: we have buffering in the submembrane space and 
//the myoplasm.
double CCell::comp_inst_buffer(double ca) {
  const double bcal = 24.0;
  const double xkcal = 7.0;
  const double srmax = 47.0;
  const double srkd = 0.6;
  const double bmem = 15.0;
  const double kmem = 0.3;
  const double bsar = 42.0;
  const double ksar = 13.0;
  double bpx = bcal*xkcal / ((xkcal + ca)*(xkcal + ca));
  double spx = srmax*srkd / ((srkd + ca)*(srkd + ca));
  double mempx = bmem*kmem / ((kmem + ca)*(kmem + ca));
  double sarpx = bsar*ksar / ((ksar + ca)*(ksar + ca));
  return 1.0 / (1.0 + bpx + spx + mempx + sarpx);
}
// --------- release-load functional dependence ----------------
double CCell::comp_Q(void) {
  double bv = (cstar - 50.) - av*cstar;
  double Qr = 0;
  if (cjsr > 50.0 && cjsr < cstar) {
    Qr = cjsr - 50.0;
  }
  else if (cjsr >= cstar) {
    Qr = av*cjsr + bv;
  }
  return cnsr*Qr / cstar;
}
double CCell::comp_dir(double po, double Qr, double rxa, double dcnsr) {
  const double ay = 0.05;
  double sparkV = exp(-ay*(v + 30)) / (1. + exp(-ay*(v + 30)));
  const double gryr = 2.58079;
  double spark_rate = gryr*po*fabs(rxa)*sparkV;
  return spark_rate*Qr - xir*(1 - taur*dcnsr / cnsr) / taur;
}
// ----------- dyadic junction dynamics ------------------------
double CCell::comp_dcp(double po, double Qr, double rxa) {
  const double gbarsr = 26841.8;// m mol/(cm C) 
  const double ax = 0.3576;
  const double gdyad = 9000.0;// m mol/(cm C) 
  double gsr = gbarsr*exp(-ax*(v + 30)) / (1.0 + exp(-ax*(v + 30)));
  double xirp = po*Qr*fabs(rxa)*gsr;

  double xicap = po*gdyad*fabs(rxa);
  const double taups = 0.5;
  return xirp + xicap - (cp - cs) / taups;
}

