#include "recorder2.h"

Recorder2::Recorder2(class Tissue *tiss, class REC2PARAM *param) {
  rcprm = param;
  X = tiss->getX();
  Y = tiss->getY();
  step = param->step;
  if (step < 1)step = 1;
  gzip = param->gzip;


  vfrm = 0;
  cfrm = 0;
#ifdef __linux__ 
  char dirtemplate[] = "/dev/shm/dirXXXXXX";
  char* ret = mkdtemp(dirtemplate);
  picdir = dirtemplate;
  picdir += "/";
#else
  picdir = "";
#endif

  savepics = false;
  picvmax = tis->cell[0].getvmax();
  picvmin = tis->cell[0].getvmin();
  piccimax = 2;
  piccimin = 0.1;

  vcolormap = COLORMAPDEFAULT;
  cicolormap = COLORMAPDEFAULT;
  tis = tiss;

  IsRectout = IsRecvout = IsRecvmovieout = IsRecekgout = false;
  if (param->filenamev != "" && param->filenamevmm != "") {
    vout = new ofstream(param->filenamev.c_str(), ios::out | ios::binary);
    vmmout = new ofstream(param->filenamevmm.c_str());
    IsRecvout = true;
  }
  if (param->filenamevmovie != "") {
    IsRecvmovieout = true;
  }
  if (param->filenamet != "") {
    tout = new ofstream(param->filenamet.c_str());
    IsRectout = true;
  }


#ifndef VOLTAGE_ONLY
  IsRecciout = IsReccnsrout = IsReccimovieout = false;
  if (param->filenamec != "" && param->filenamecmm != "") {
    ciout = new ofstream(param->filenamec.c_str(), ios::out | ios::binary);
    cimmout = new ofstream(param->filenamecmm.c_str());
    IsRecciout = true;
  }
  if (param->filenamecmovie != "") {
    IsReccimovieout = true;
  }
  if (param->filenamecnsr != "" && param->filenamecnsrmm != "") {
    cnsrout = new ofstream(param->filenamecnsr.c_str(), ios::out | ios::binary);
    cnsrmmout = new ofstream(param->filenamecnsrmm.c_str());
    IsReccnsrout = true;
  }
#endif

  if (param->filenameekg != "") {
    ekgout = new ofstream(param->filenameekg.c_str());
    IsRecekgout = true;
  }
}

Recorder2::~Recorder2(void) {
  if (IsRecvout) {
    vout->close(); vmmout->close(); delete vout; delete vmmout;
    if (gzip){ char temp[filenamemaxchar]; sprintf(temp, "xz -9f %s", rcprm->filenamev.c_str()); int res = system(temp); if (res)cout << res << endl; }
  }
  if (IsRecvmovieout) {
    char temp[2000];
    if (rcprm->pictureformat == 3)
      sprintf(temp, "ffmpeg -i %svframe%%08d.jpg -preset placebo -y -r 30 %s %s;", picdir.c_str(), rcprm->ffmpegoptions.c_str(), rcprm->filenamevmovie.c_str());
    else if (rcprm->pictureformat == 2)
      sprintf(temp, "ffmpeg -i %svframe%%08d.webp -preset placebo -y -r 30 %s %s;", picdir.c_str(), rcprm->ffmpegoptions.c_str(), rcprm->filenamevmovie.c_str());
    else if (rcprm->pictureformat == 1)
      sprintf(temp, "ffmpeg -i %svframe%%08d.png -preset placebo -y -r 30 %s %s;", picdir.c_str(), rcprm->ffmpegoptions.c_str(), rcprm->filenamevmovie.c_str());
    else
      sprintf(temp, "ffmpeg -i %svframe%%08d.ppm -preset placebo -y -r 30 %s %s;", picdir.c_str(), rcprm->ffmpegoptions.c_str(), rcprm->filenamevmovie.c_str());
    int res = system(temp);
    if (res)cout << res << endl;
    if (savepics == false) {
      sprintf(temp, "rm %svframe*.*", picdir.c_str());
      res = system(temp);
      if (res)cout << res << endl;
    }
  }
  if (IsRectout){ tout->close(); delete tout; }
  if (IsRecekgout) {
    ekgout->close(); delete ekgout;
    if (gzip){ char temp[filenamemaxchar]; sprintf(temp, "xz -9f %s", rcprm->filenameekg.c_str()); int res = system(temp); if (res)cout << res << endl; }
  }

#ifndef VOLTAGE_ONLY
  if (IsRecciout) {
    ciout->close(); cimmout->close(); delete ciout; delete cimmout;
    if (gzip){ char temp[filenamemaxchar]; sprintf(temp, "xz -9f %s", rcprm->filenamec.c_str()); int res = system(temp); if (res)cout << res << endl; }
  }
  if (IsReccimovieout) {
    char temp[2000];
    if (rcprm->pictureformat == 3)
      sprintf(temp, "ffmpeg -i %scframe%%08d.jpg -preset placebo -y -r 30 %s %s;", picdir.c_str(), rcprm->ffmpegoptions.c_str(), rcprm->filenamecmovie.c_str());
    else if (rcprm->pictureformat == 2)
      sprintf(temp, "ffmpeg -i %scframe%%08d.webp -preset placebo -y -r 30 %s %s;", picdir.c_str(), rcprm->ffmpegoptions.c_str(), rcprm->filenamecmovie.c_str());
    else if (rcprm->pictureformat == 1)
      sprintf(temp, "ffmpeg -i %scframe%%08d.png -preset placebo -y -r 30 %s %s;", picdir.c_str(), rcprm->ffmpegoptions.c_str(), rcprm->filenamecmovie.c_str());
    else
      sprintf(temp, "ffmpeg -i %scframe%%08d.ppm -preset placebo -y -r 30 %s %s;", picdir.c_str(), rcprm->ffmpegoptions.c_str(), rcprm->filenamecmovie.c_str());
    int res = system(temp);
    if (res)cout << res << endl;
    if (savepics == false) {
      sprintf(temp, "rm %scframe*.*", picdir.c_str());
      res = system(temp);
      if (res)cout << res << endl;
    }
  }
  if (IsReccnsrout) {
    cnsrout->close(); cnsrmmout->close(); delete cnsrout; delete cnsrmmout;
    if (gzip){ char temp[filenamemaxchar]; sprintf(temp, "xz -9f %s", rcprm->filenamecnsr.c_str()); int res = system(temp); if (res)cout << res << endl; }
  }
#endif
}

void Recorder2::recdata(double t) {
  if (IsRectout)*tout << t << endl;

  if (IsRecvout) {
    double vmax = -FLT_MAX;
    double vmin = FLT_MAX;
    for (int i = 0; i < X; i += step) {
      for (int j = 0; j<Y; j += step) {
        if (tis->cell[i*Y + j].v>vmax)vmax = tis->cell[i*Y + j].v;
        if (tis->cell[i*Y + j].v < vmin)vmin = tis->cell[i*Y + j].v;
      }
    }
    if (vmax - vmin < 1)vmax = vmin + 1;
    double vdiff = vmax - vmin;


    *vmmout << vmax << "\t" << vmin << endl;
    for (int i = 0; i < X; i += step) {
      for (int j = 0; j < Y; j += step) {
        unsigned char urectmp;
        urectmp = (unsigned char)((tis->cell[i*Y + j].v - vmin)*(UCHAR_MAX / vdiff));
        vout->write((char*)&urectmp, sizeof(unsigned char));
      }
    }
    vout->flush();
  }


  if (IsRecvmovieout) {
    char fname[255];
    sprintf(fname, "%svframe%08d.ppm", picdir.c_str(), vfrm);
    ofstream osfrm(fname, ofstream::binary);

    osfrm << "P6\n";
    osfrm << "# vframe\n";
    osfrm << X << " " << Y << "\n";
    osfrm << "255\n";
    char* tmp = new char[X*Y * 3];

    if (vcolormap == COLORMAPJET) {
      double picvdiff = picvmax - picvmin;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < X; i++) {
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector always
#endif
        for (int j = 0; j < Y; j++) {
          double normal = (tis->cell[i*Y + j].v - picvmin) / picvdiff;
          if (normal < 0)normal = 0;
          if (normal > 1)normal = 1;

          char r = 0;
          char g = 0;
          char b = 0;

          if (normal > 0.375) {
            if (normal<0.625)
              r = (char)((normal*4.0 - 1.5) * 255);
            else if (normal>0.875)
              r = (char)((4.5 - normal*4.0) * 255);
            else
              r = (char)(255);
          }

          if (normal > 0.125) {
            if (normal < 0.375)
              g = (char)((4.0*normal - 0.5) * 255);
            else if (normal < 0.625)
              g = (char)(255);
            else if (normal < 0.875)
              g = (char)((3.5 - normal*4.0) * 255);
          }

          if (normal < 0.125)
            b = (char)((0.5 + normal*4.0) * 255);
          else if (normal < 0.375)
            b = (char)(255);
          else if (normal < 0.625)
            b = (char)((2.5 - normal*4.0) * 255);

          int loc = (i*Y + j) * 3;
          tmp[loc] = r;
          tmp[loc + 1] = g;
          tmp[loc + 2] = b;
        }
      }
      osfrm.write((char*)tmp, sizeof(char)*X*Y * 3);
    }
    else if (vcolormap == COLORMAPHOT) {
      double picvdiff = picvmax - picvmin;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < X; i++) {
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector always
#endif
        for (int j = 0; j < Y; j++) {
          double normal = (tis->cell[i*Y + j].v - picvmin) / picvdiff;
          if (normal < 0)normal = 0;
          if (normal > 1)normal = 1;

          char r = (char)(255);
          char g = 0;
          char b = 0;

          if (normal < 0.375) {
            r = (char)((normal*1.0 / 0.375) * 255);
          }
          else if (normal < 0.750) {
            g = (char)((normal*1.0 / 0.375 - 1.0) * 255);
          }
          else {
            b = (char)((normal*4.0 - 3.0) * 255);
            g = (char)(255);
          }

          int loc = (i*Y + j) * 3;
          tmp[loc] = r;
          tmp[loc + 1] = g;
          tmp[loc + 2] = b;
        }
      }
      osfrm.write((char*)tmp, sizeof(char)*X*Y * 3);
    }
    else if (vcolormap == COLORMAPGRAY) {
      double picvdiff = picvmax - picvmin;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < X; i++) {
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector always
#endif
        for (int j = 0; j < Y; j++) {
          double normal = (tis->cell[i*Y + j].v - picvmin) / picvdiff;
          if (normal < 0)normal = 0;
          if (normal > 1)normal = 1;

          char r = (char)(normal * 255);
          char g = (char)(normal * 255);
          char b = (char)(normal * 255);

          int loc = (i*Y + j) * 3;
          tmp[loc] = r;
          tmp[loc + 1] = g;
          tmp[loc + 2] = b;
        }
      }
      osfrm.write((char*)tmp, sizeof(char)*X*Y * 3);
    }
    else {
      double picvdiff = picvmax - picvmin;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < X; i++) {
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector always
#endif
        for (int j = 0; j < Y; j++) {
          double normal = (tis->cell[i*Y + j].v - picvmin) / picvdiff;

          char r = 0;
          char g = (char)(255);
          char b = 0;
          if (normal > 1) {
            r = (char)(255);
            g = 0;
          }
          else if (normal > 0.75) {
            r = (char)(255);
            g = (char)((-4 * normal + 4) * 255);
          }
          else if (normal > 0.5) {
            r = (char)((4 * normal - 2) * 255);
          }
          else if (normal > 0.25) {
            b = (char)((-4 * normal + 2) * 255);
          }
          else if (normal > 0) {
            g = (char)((4 * normal) * 255);
            b = (char)(255);
          }
          else {
            g = 0;
            b = (char)(255);
          }
          int loc = (i*Y + j) * 3;
          tmp[loc] = r;
          tmp[loc + 1] = g;
          tmp[loc + 2] = b;
        }
      }
      osfrm.write((char*)tmp, sizeof(char)*X*Y * 3);
    }

    delete[] tmp;

    if (rcprm->pictureformat > 0)
    {
      char cmdtmp[2550];
      if (rcprm->pictureformat == 1)
        sprintf(cmdtmp, "convert %svframe%08d.ppm -depth 24 %svframe%08d.png;rm %svframe%08d.ppm &", picdir.c_str(), vfrm, picdir.c_str(), vfrm, picdir.c_str(), vfrm);
      if (rcprm->pictureformat == 2)
        sprintf(cmdtmp, "convert %svframe%08d.ppm -depth 24 %svframe%08d.png;cwebp %svframe%08d.png -lossless -quiet -o %svframe%08d.webp;rm %svframe%08d.ppm;rm %svframe%08d.png &", picdir.c_str(), vfrm, picdir.c_str(), vfrm, picdir.c_str(), vfrm, picdir.c_str(), vfrm, picdir.c_str(), vfrm, picdir.c_str(), vfrm);
      if (rcprm->pictureformat == 3)
        sprintf(cmdtmp, "convert %svframe%08d.ppm -depth 24 %svframe%08d.jpg;rm %svframe%08d.ppm &", picdir.c_str(), vfrm, picdir.c_str(), vfrm, picdir.c_str(), vfrm);
      int res = system(cmdtmp);
      if (res)cout << res << endl;
    }
    vfrm++;
  }


#ifndef VOLTAGE_ONLY
  if (IsRecciout) {
    double caimax = -FLT_MAX;
    double caimin = FLT_MAX;
    for (int i = 0; i < X; i += step) {
      for (int j = 0; j<Y; j += step) {
        if (tis->cell[i*Y + j].ci>caimax)caimax = tis->cell[i*Y + j].ci;
        if (tis->cell[i*Y + j].ci < caimin)caimin = tis->cell[i*Y + j].ci;
      }
    }
    if (caimax - caimin < 0.001)caimax = caimin + 0.001;
    double cidiff = caimax - caimin;
    *cimmout << caimax << "\t" << caimin << endl;
    for (int i = 0; i < X; i += step) {
      for (int j = 0; j < Y; j += step) {
        unsigned char urectmp;
        urectmp = (unsigned char)((tis->cell[i*Y + j].ci - caimin)*(UCHAR_MAX / cidiff));
        ciout->write((char *)&urectmp, sizeof(unsigned char));
      }
    }
    ciout->flush();
  }


  if (IsReccimovieout) {
    char fname[255];
    sprintf(fname, "%scframe%08d.ppm", picdir.c_str(), cfrm);
    ofstream osfrm(fname, ofstream::binary);

    osfrm << "P6\n";
    osfrm << "# cframe\n";
    osfrm << X << " " << Y << "\n";
    osfrm << "255\n";
    char* tmp = new char[X*Y * 3];

    if (cicolormap == COLORMAPJET) {
      double piccidiff = piccimax - piccimin;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < X; i++) {
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector always
#endif
        for (int j = 0; j < Y; j++) {
          double normal = (tis->cell[i*Y + j].ci - piccimin) / piccidiff;
          if (normal < 0)normal = 0;
          if (normal > 1)normal = 1;

          char r = 0;
          char g = 0;
          char b = 0;

          if (normal > 0.375) {
            if (normal<0.625)
              r = (char)((normal*4.0 - 1.5) * 255);
            else if (normal>0.875)
              r = (char)((4.5 - normal*4.0) * 255);
            else
              r = (char)(255);
          }

          if (normal > 0.125) {
            if (normal < 0.375)
              g = (char)((4.0*normal - 0.5) * 255);
            else if (normal < 0.625)
              g = (char)(255);
            else if (normal < 0.875)
              g = (char)((3.5 - normal*4.0) * 255);
          }

          if (normal < 0.125)
            b = (char)((0.5 + normal*4.0) * 255);
          else if (normal < 0.375)
            b = (char)(255);
          else if (normal < 0.625)
            b = (char)((2.5 - normal*4.0) * 255);

          int loc = (i*Y + j) * 3;
          tmp[loc] = r;
          tmp[loc + 1] = g;
          tmp[loc + 2] = b;
        }
      }
      osfrm.write((char*)tmp, sizeof(char)*X*Y * 3);
    }
    else if (cicolormap == COLORMAPHOT) {
      double piccidiff = piccimax - piccimin;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < X; i++) {
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector always
#endif
        for (int j = 0; j < Y; j++) {
          double normal = (tis->cell[i*Y + j].ci - piccimin) / piccidiff;
          if (normal < 0)normal = 0;
          if (normal > 1)normal = 1;

          char r = (char)(255);
          char g = 0;
          char b = 0;

          if (normal < 0.375) {
            r = (char)((normal*1.0 / 0.375) * 255);
          }
          else if (normal < 0.750) {
            g = (char)((normal*1.0 / 0.375 - 1.0) * 255);
          }
          else {
            b = (char)((normal*4.0 - 3.0) * 255);
            g = (char)(255);
          }

          int loc = (i*Y + j) * 3;
          tmp[loc] = r;
          tmp[loc + 1] = g;
          tmp[loc + 2] = b;
        }
      }
      osfrm.write((char*)tmp, sizeof(char)*X*Y * 3);
    }
    else if (cicolormap == COLORMAPGRAY) {
      double piccidiff = piccimax - piccimin;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < X; i++) {
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector always
#endif
        for (int j = 0; j < Y; j++) {
          double normal = (tis->cell[i*Y + j].ci - piccimin) / piccidiff;
          if (normal < 0)normal = 0;
          if (normal > 1)normal = 1;

          char r = (char)(normal * 255);
          char g = (char)(normal * 255);
          char b = (char)(normal * 255);

          int loc = (i*Y + j) * 3;
          tmp[loc] = r;
          tmp[loc + 1] = g;
          tmp[loc + 2] = b;
        }
      }
      osfrm.write((char*)tmp, sizeof(char)*X*Y * 3);
    }
    else {
      double piccidiff = piccimax - piccimin;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < X; i++) {
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector always
#endif
        for (int j = 0; j < Y; j++) {
          double normal = (tis->cell[i*Y + j].ci - piccimin) / piccidiff;

          char r = 0;
          char g = (char)(255);
          char b = 0;
          if (normal > 1) {
            r = (char)(255);
            g = 0;
          }
          else if (normal > 0.75) {
            r = (char)(255);
            g = (char)((-4 * normal + 4) * 255);
          }
          else if (normal > 0.5) {
            r = (char)((4 * normal - 2) * 255);
          }
          else if (normal > 0.25) {
            b = (char)((-4 * normal + 2) * 255);
          }
          else if (normal > 0) {
            g = (char)((4 * normal) * 255);
            b = (char)(255);
          }
          else {
            g = 0;
            b = (char)(255);
          }
          int loc = (i*Y + j) * 3;
          tmp[loc] = r;
          tmp[loc + 1] = g;
          tmp[loc + 2] = b;
        }
      }
      osfrm.write((char*)tmp, sizeof(char)*X*Y * 3);
    }

    delete[] tmp;

    if (rcprm->pictureformat > 0) {
      char cmdtmp[2550];
      if (rcprm->pictureformat == 1)
        sprintf(cmdtmp, "convert %scframe%08d.ppm -depth 24 %scframe%08d.png;rm %scframe%08d.ppm &", picdir.c_str(), cfrm, picdir.c_str(), cfrm, picdir.c_str(), cfrm);
      if (rcprm->pictureformat == 2)
        sprintf(cmdtmp, "convert %scframe%08d.ppm -depth 24 %scframe%08d.png;cwebp %scframe%08d.png -lossless -quiet -o %scframe%08d.webp;rm %scframe%08d.ppm;rm %scframe%08d.png &", picdir.c_str(), cfrm, picdir.c_str(), cfrm, picdir.c_str(), cfrm, picdir.c_str(), cfrm, picdir.c_str(), cfrm, picdir.c_str(), cfrm);
      if (rcprm->pictureformat == 3)
        sprintf(cmdtmp, "convert %scframe%08d.ppm -depth 24 %scframe%08d.jpg;rm %scframe%08d.ppm &", picdir.c_str(), cfrm, picdir.c_str(), cfrm, picdir.c_str(), cfrm);
      int res = system(cmdtmp);
      if (res)cout << res << endl;
    }
    cfrm++;
  }



  if (IsReccnsrout) {
    double cansrmax = -FLT_MAX;
    double cansrmin = FLT_MAX;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < X; i += step) {
      for (int j = 0; j<Y; j += step) {
        if (tis->cell[i*Y + j].cnsr>cansrmax)cansrmax = tis->cell[i*Y + j].cnsr;
        if (tis->cell[i*Y + j].cnsr < cansrmin)cansrmin = tis->cell[i*Y + j].cnsr;
      }
    }
    if (cansrmax - cansrmin < 0.001)cansrmax = cansrmin + 0.001;
    double cidiff = cansrmax - cansrmin;
    *cnsrmmout << cansrmax << "\t" << cansrmin << endl;
    for (int i = 0; i < X; i += step) {
      for (int j = 0; j < Y; j += step) {
        unsigned char urectmp;
        urectmp = (unsigned char)((tis->cell[i*Y + j].cnsr - cansrmin)*(UCHAR_MAX / cidiff));
        cnsrout->write((char *)&urectmp, sizeof(unsigned char));
      }
    }
    cnsrout->flush();
  }
#endif
}



void Recorder2::recekg(double t, double posx, double posy, double posz, bool first, bool last) {

  double dx = tis->getdx();
  double ekg = 0;
  double Dfudx = tis->getdif() / (dx*dx);
  for (int i = 1; i < X - 1; i++) {
    for (int j = 1; j < Y - 1; j++) {
      double ix = i*dx;
      double iy = j*dx;
      double r = sqrt((ix - posx)*(ix - posx) + (iy - posy)*(iy - posy) + posz*posz);

      double t_ekgx = (tis->cell[(i + 1)*Y + j].v - tis->cell[(i - 1)*Y + j].v)*Dfudx*(posx - ix);
      double t_ekgy = (tis->cell[i*Y + j + 1].v - tis->cell[i*Y + j - 1].v)*Dfudx*(posy - iy);

      ekg += (t_ekgx + t_ekgy) / (r*r*r);
    }
  }
  if (IsRecekgout) {
    if (first)*ekgout << t;
    *ekgout << "\t" << ekg;
    if (last)*ekgout << endl;
  }
}


