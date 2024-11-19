/***************************************************************************
 *
 *   Copyright (C) 2012 by Cees Bassa
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <cpgplot.h>

#include <iostream>
#include <fstream>

#include "Pulsar/psrchive.h"

#include "Pulsar/Archive.h"
#include "Pulsar/Profile.h"
#include "Pulsar/Integration.h"

#include "Pulsar/RemoveVariableBaseline.h"

#include "strutil.h"

using namespace std;
using namespace Pulsar;

enum cursor_type {
  horizontal,
  vertical,
  both
};

enum zap_type {
  freq_zap,
  time_zap,
  both_zap
};

enum plot_type {
  freq_plot,
  time_plot,
  both_plot,
  prof_plot,
  wt_plot
};

struct zap_range {
  int chan0;
  int chan1;
  int sub0;
  int sub1;
  enum zap_type type;
};

struct image {
  int nx,ny;
  float *z;
  float tr[6];
};

void swap_zap_range(struct zap_range *z) 
{
  int tmp;
  if (z->chan0 > z->chan1) { tmp=z->chan0; z->chan0=z->chan1; z->chan1=tmp; }
  if (z->sub0 > z->sub1) { tmp=z->sub0; z->sub0=z->sub1; z->sub1=tmp; }
}

void apply_zap_range(Archive *arch, struct zap_range *z,bool undo=false, const Archive *orig=NULL) {
  swap_zap_range(z);
  int nsub = arch->get_nsubint();
  int nchan = arch->get_nchan();

  // Keep in range
  if (z->sub0 < 0) z->sub0 = 0;
  if (z->sub1 >= nsub) z->sub1 = nsub-1;
  if (z->chan0 < 0) z->chan0 = 0;
  if (z->chan1 >= nchan) z->chan1 = nchan-1;

  for (int isub=z->sub0; isub<=z->sub1; isub++) {
    Reference::To<Integration> subint = arch->get_Integration(isub);
    for (int ichan=z->chan0; ichan<=z->chan1; ichan++) {
      if (undo && orig!=NULL)
        subint->set_weight(ichan, 
            orig->get_Integration(isub)->get_weight(ichan));
      else
        subint->set_weight(ichan,0.0);
    }
  }
}

void zap_freq_range(Archive *arch, struct zap_range *z) 
{
  z->sub0 = 0;
  z->sub1 = arch->get_nsubint()-1;
  apply_zap_range(arch, z);
}

void zap_subint_range(Archive *arch, struct zap_range *z) {
  swap_zap_range(z);
  int nsub = arch->get_nsubint();
  int nchan = arch->get_nchan();
  z->chan0 = 0;
  z->chan1 = nchan;
  if (z->sub0 < 0) z->sub0 = 0;
  if (z->sub1 >= nsub) z->sub1 = nsub-1;
  for (int isub=z->sub0; isub<=z->sub1; isub++) 
    arch->get_Integration(isub)->uniform_weight(0.0);
}

void apply_zap(Archive *arch, struct zap_range *z) {
  if (z->type==freq_zap) 
    zap_freq_range(arch, z);
  else if (z->type==time_zap) 
    zap_subint_range(arch, z);
  else if (z->type==both_zap) 
    apply_zap_range(arch, z);
}

// Get extrema
void minmax(float x[],int n,float *xmin,float *xmax)
{
  int i,flag=0;
  float c;

  for (i=0;i<n;i++) {
    if (flag==0) {
      *xmin= *xmax=x[i];
      flag=1;
    }
    if (x[i]< *xmin) *xmin=x[i];
    if (x[i]> *xmax) *xmax=x[i];
  }
  c=0.1*( *xmax- *xmin);

  *xmin-=c;
  *xmax+=c;

  return;
}

// Create yfp image
struct image create_yfp_image(Pulsar::Archive* archive)
{
  int isub,ibin,ipol,ichan,k;
  int nsub,nbin,npol,nchan;
  struct image img;
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};

  // Get image sizes
  nsub = archive->get_nsubint();
  nbin = archive->get_nbin();
  nchan = archive->get_nchan();
  npol = archive->get_npol();

  // Allocate
  img.nx=nbin;
  img.ny=nsub;
  img.z=(float *) malloc(sizeof(float)*nbin*nsub);

  // Set transformation matrix
  for (k=0;k<6;k++)
    img.tr[k]=tr[k];
  img.tr[0]=-0.5/(float) nbin;
  img.tr[1]=1.0/(float) nbin;
  
  // Set to zero
  for (k=0;k<nbin*nsub;k++)
    img.z[k]=0.0;

  // Loop over archive
  for (isub=0;isub<nsub;isub++) {
    for (ichan=0;ichan<nchan;ichan++) {
      for (ipol=0;ipol<npol;ipol++) {
	Pulsar::Profile* profile = archive->get_Profile (isub,ipol,ichan);
	float* data=profile->get_amps();
	float wt=profile->get_weight();
	
	for (ibin=0;ibin<nbin;ibin++) {
	  k=nbin*isub+ibin;
	  if (wt!=0.0) 
	    img.z[k]+=data[ibin]*wt;
	}
      }
    }
  }

  return img;
}

// Create gtp image
struct image create_gtp_image(Pulsar::Archive* archive)
{
  int isub,ibin,ipol,ichan,k;
  int nsub,nbin,npol,nchan;
  struct image img;
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};

  // Get image sizes
  nsub = archive->get_nsubint();
  nbin = archive->get_nbin();
  nchan = archive->get_nchan();
  npol = archive->get_npol();

  // Allocate
  img.nx=nbin;
  img.ny=nchan;
  img.z=(float *) malloc(sizeof(float)*nbin*nchan);

  // Set transformation matrix
  for (k=0;k<6;k++)
    img.tr[k]=tr[k];
  img.tr[0]=-0.5/(float) nbin;
  img.tr[1]=1.0/(float) nbin;

  for (k=0;k<nbin*nchan;k++)
    img.z[k]=0.0;
  
  // Loop over archive
  for (ichan=0;ichan<nchan;ichan++) {
    for (isub=0;isub<nsub;isub++) {
      for (ipol=0;ipol<npol;ipol++) {
	Pulsar::Profile* profile = archive->get_Profile (isub,ipol,ichan);
	float* data=profile->get_amps();
	float wt=profile->get_weight();
	
	for (ibin=0;ibin<nbin;ibin++) {
	  k=nbin*ichan+ibin;
	  if (wt!=0.0) 
	    img.z[k]+=data[ibin]*wt;
	}
      }
    }
  }

  return img;
}

// Create dynamic spectrum
struct image create_ds_image(Pulsar::Archive* archive)
{
  int isub,ichan,ibin,k;
  int nsub,nchan,nbin;
  struct image img;
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};

  // Get image sizes
  nsub = archive->get_nsubint();
  nchan = archive->get_nchan();
  nbin = archive->get_nbin();

  // Allocate
  img.nx=nsub;
  img.ny=nchan;
  img.z=(float *) malloc(sizeof(float)*nsub*nchan);

  // Set transformation matrix
  for (k=0;k<6;k++)
    img.tr[k]=tr[k];
  
  // Loop over archive
  for (ichan=0; ichan < nchan; ichan++) {
    for (isub=0; isub < nsub; isub++) {
      Pulsar::Profile* profile = archive->get_Profile (isub,0,ichan);
      float* data=profile->get_amps();
      float wt=profile->get_weight();
      k=nsub*ichan+isub;

      if (wt!=0.0) {
	float s1=0.0;
	float s2=0.0;
	// compute variance and mean
	for (ibin=0;ibin<nbin;ibin++) {
	  s1+=data[ibin];
	  s2+=data[ibin]*data[ibin];
	}
	float mean=s1/(float) nbin;
	float var=(s2-s1*mean)/(float) (nbin-1);
	img.z[k]=var;
      } else {
	img.z[k]=0.0;
      }
    }
  }

  return img;
}

// Create weight image
struct image create_wt_image(Pulsar::Archive* archive)
{
  int isub,ichan,ibin,k;
  int nsub,nchan,nbin;
  struct image img;
  float tr[]={-0.5,1.0,0.0,-0.5,0.0,1.0};

  // Get image sizes
  nsub = archive->get_nsubint();
  nchan = archive->get_nchan();
  nbin = archive->get_nbin();

  // Allocate
  img.nx=nsub;
  img.ny=nchan;
  img.z=(float *) malloc(sizeof(float)*nsub*nchan);

  // Set transformation matrix
  for (k=0;k<6;k++)
    img.tr[k]=tr[k];
  
  // Loop over archive
  for (ichan=0; ichan < nchan; ichan++) {
    for (isub=0; isub < nsub; isub++) {
      Pulsar::Profile* profile = archive->get_Profile (isub,0,ichan);
      float wt=profile->get_weight();
      k=nsub*ichan+isub;

      img.z[k]=wt;
    }
  }

  return img;
}

// Plot a fully scrunched profile
void plot_profile(struct image img)
{
  int i,j,k;
  float *z,zmin,zmax;
  
  // Allocate
  z=(float *) malloc(sizeof(float)*img.nx);
  
  // Set to zero
  for (k=0;k<img.nx;k++)
    z[k]=0.0;
  
  // Loop over image
  for (i=0;i<img.nx;i++) {
    for (j=0;j<img.ny;j++) {
      k=img.nx*j+i;
      z[i]+=img.z[k];
    }
  }

  // Scale back
  for (k=0;k<img.nx;k++)
    z[k]/=(float) img.nx;
  
  minmax(z,img.nx,&zmin,&zmax);
  
  cpgswin(0.0,1.0,zmin,zmax);
  cpgbox("BCTSNI",0.,0,"BCTSNI",0.,0);
  cpglab("Pulse Phase","Intensity"," ");

  for (i=0;i<img.nx;i++) {
    if (i==0)
      cpgmove((float) i/(float) img.nx,z[i]);
    else
      cpgdraw((float) i/(float) img.nx,z[i]);
  }

  free(z);

  return;
}

// Plot an image
void plot_image(struct image img,bool dynamic_levels)
{
  int i,j,k,flag=0;
  float zmin,zmax;
  float xmin,xmax,ymin,ymax;
  int imin,imax,jmin,jmax;
  float s0=0.0,s1=0.0,mean,std;
  int sn=0;

  // Get plot limits
  cpgqwin(&xmin,&xmax,&ymin,&ymax);
  imin=(int) ((xmin-img.tr[0])/img.tr[1]-0.5);
  imax=(int) ((xmax-img.tr[0])/img.tr[1]+0.5);
  jmin=(int) ((ymin-img.tr[3])/img.tr[5]-0.5);
  jmax=(int) ((ymax-img.tr[3])/img.tr[5]+0.5);

  // Keep in range
  if (imin<0) imin=0;
  if (imax>img.nx-1) imax=img.nx-1;
  if (jmin<0) jmin=0;
  if (jmax>img.ny-1) jmax=img.ny-1;

  // Determine extrema
  for (i=imin;i<imax;i++) {
    for (j=jmin;j<jmax;j++) {
      k=img.nx*j+i;
      if (flag==0) {
	zmin=img.z[k];
	zmax=img.z[k];
	flag=1;
      }
      if (img.z[k]<zmin) zmin=img.z[k];
      if (img.z[k]>zmax) zmax=img.z[k];
      s0+=img.z[k];
      s1+=img.z[k]*img.z[k];
      sn++;
    }
  }
  if (dynamic_levels==true) {
    mean=s0/(float) sn;
    std=sqrt((s1-s0*mean)/(float) (sn-1));
    
    //    zmin=mean-4.0*std;
    zmax=mean+6.0*std;
  }

  // Plot image
  cpgimag(img.z,img.nx,img.ny,1,img.nx,1,img.ny,zmin,zmax,img.tr);

  return;
}


#define PROG "psrzap"

/* Extension for output file */
string output_ext = "zap";

void usage() {
  cout << "psrzap - Interactive RFI zapper using off-pulse dynamic spectra" 
    << endl
    << "Usage:  psrzap (filename)" << endl
    << "Command line options:" << endl
    << "  -h     Show this help screen" << endl
    << "  -v     Verbose mode" << endl
    << "  -e ext Extension for output file (default: " << output_ext 
      << ")" << endl
    << endl
    << "  Press 'h' during program for interactive usage info:" << endl
    << endl;
}

/* Interactive commands */
#define CMD_QUIT 'q'
#define CMD_HELP 'h'
#define CMD_SAVE 's'
//#define CMD_PRINT 'P'
#define CMD_FREQMODE 'f'
#define CMD_TIMEMODE 't'
#define CMD_BOTHMODE 'b'
#define CMD_BOTH '3'
#define CMD_TIME '1'
#define CMD_FREQ '2'
#define CMD_PROF '4'
#define CMD_WT '5'
//#define CMD_POL 'p'
//#define CMD_VAR 'v'
//#define CMD_LOG 'l'
#define CMD_UNDO 'u'
#define CMD_UNZOOM 'r'
#define CMD_UPDATE 'D'
#define CMD_DEDISPERSE 'd'
#define CMD_PSRSH 'w'
#define CMD_LEVELS 'e'
//#define CMD_BASELINE 'B'
//#define CMD_METHOD 'm'
void usage_interactive() {
  cout << "pzrzap interactive commands" << endl
       << endl
       << "Mouse:" << endl
       << "  Left-click selects the start of a range" << endl
       << "    then left-click again to zoom, or right-click to zap." << endl
       << "  Right-click zaps current cursor location." << endl
       << "  Middle-click (or 'D') updates the diagnostic plots." << endl
       << "Keyboard:" << endl
       << "  " << CMD_HELP     << "  Show this help" << endl
       << "  " << CMD_FREQMODE << "  Use frequency select mode" << endl
       << "  " << CMD_TIMEMODE << "  Use time select mode" << endl
       << "  " << CMD_BOTHMODE << "  Use time/freq select mode" << endl
       << "  " << CMD_TIME << "  Show pulse phase vs sub-integration" << endl
       << "  " << CMD_WT << "  Show profile weights" << endl
       << "  " << CMD_FREQ << "  Show pulse phase vs channel" << endl
       << "  " << CMD_BOTH << "  Show dynamic spectrum" << endl
       << "  " << CMD_PROF << "  Show profile" << endl
       << "  " << CMD_DEDISPERSE << "  Toggle dedispersion" <<endl
       << "  " << CMD_UPDATE   << "  Update diagnostic plots" << endl
    //    << "  " << CMD_VAR      << "  Switch plot between variance and mean" << endl
    //    << "  " << CMD_LOG      << "  Toggle logscale" << endl
    //    << "  " << CMD_METHOD   << "  Toggle method" << endl
    //    << "  " << CMD_POL      << "  Switch to next polarization" << endl
    //    << "  " << CMD_BASELINE << "  Toggle variable baseline removal" << endl
       << "  " << CMD_UNDO     << "  Undo last zap command" << endl 
       << "  " << CMD_UNZOOM   << "  Reset zoom and update plots" << endl
       << "  " << CMD_SAVE     << "  Save zapped version as (filename)." 
       << output_ext << endl
    //    << "  " << CMD_PRINT    << "  Print equivalent paz command" << endl
       << "  " << CMD_PSRSH    << "  Generate PSRSH script" << endl
       << "  " << CMD_QUIT     << "  Exit program" << endl
       << endl;
}

int main(int argc, char *argv[]) {

  /* Process any args */
  int opt=0;
  int verb=0;

  string expression;

  while ((opt=getopt(argc,argv,"hve:E:"))!=-1) {
    switch (opt) {

      case 'v':
        verb++;
        Archive::set_verbosity(verb);
        break;
        
      case 'e':
        output_ext.assign(optarg);
        break;

      case 'E':
	expression = optarg;
	break;

      case 'h':
      default:
        usage();
        usage_interactive();
        exit(0);
        break;

    }
  }

  if (optind==argc) {
    usage();
    cerr << PROG ": No filename given" << endl;
    exit(-1);
  }

  // Load file 
  string filename = argv[optind];
  Reference::To<Archive> orig_arch = Archive::load(filename);
  Reference::To<Archive> arch = orig_arch->clone();
  arch->dedisperse();
  arch->remove_baseline();
  string output_filename = replace_extension(filename, output_ext);
  string psrsh_filename = replace_extension(filename,"psh");

  // Create profile plots
  struct image fp,tp,ds,wt;

  // Color table
  float hl[]={0.0,0.2,0.4,0.6,1.0},hr[]={0.0,0.5,1.0,1.0,1.0},hg[]={0.0,0.0,0.5,1.0,1.0},hb[]={0.0,0.0,0.0,0.3,1.0};

  // Image coordinates
  float x0,y0,x1,y1;
  float fmin,fmax,tmin,tmax;
  char ch='\0';
  int mode=0;
  float freqmin,freqmax,freq,bw;

  // Boolean flags
  bool redraw=true,click=false,replot=true,unzoom=true,dynamic_levels=false;
  enum cursor_type curs=horizontal;
  enum plot_type plot=time_plot;

  // Set up zap structs
  struct zap_range zap;
  vector<struct zap_range> zap_list;
  int izap=0,ipol=0;

  freq=arch->get_centre_frequency();
  bw=arch->get_bandwidth();

  // Set up plot
  cpgopen("/xw");
  cpgask(0);
  cpgsfs(2);
  cpgctab(hl,hr,hg,hb,5,1.0,0.5);

  do {
    // Redraw images
    if (redraw) {
      fp=create_yfp_image(arch);
      tp=create_gtp_image(arch);
      ds=create_ds_image(arch);
      wt=create_wt_image(arch);
    }

    // Get image limits
    if (unzoom) {
      tmin=0.0;
      tmax=(float) ds.nx;
      fmin=0.0;
      fmax=(float) ds.ny;
      unzoom=false;
    }
    freqmin=freq-0.5*bw;
    freqmax=freq+0.5*bw;
    freqmin=freq+(fmin/ds.ny-0.5)*bw;
    freqmax=freq+(fmax/ds.ny-0.5)*bw;

    // Plot images
    if (replot==true) {
      if (plot==time_plot) {
	cpgeras();
	cpgswin(0.0,1.0,tmin,tmax);
	cpglab("Pulse phase","Sub Integration","");
	plot_image(fp,dynamic_levels);
	cpgbox("BCTSNI",0.,0,"BCTSNI",0.,0);
      } else if (plot==freq_plot) {
	cpgeras();
	cpgswin(0.0,1.0,freqmin,freqmax);
	cpgbox("BCTSNI",0.,0,"CTSMI",0.,0);
	cpgswin(0.0,1.0,fmin,fmax);
	cpglab("Pulse phase","Channel","");
	cpgmtxt("R",3.0,0.5,0.5,"Frequency (MHz)");
	plot_image(tp,dynamic_levels);
	cpgbox("BCTSNI",0.,0,"BTSNI",0.,0);
      } else if (plot==both_plot) {
	cpgeras();
	cpgswin(tmin,tmax,freqmin,freqmax);
	cpgbox("BCTSNI",0.,0,"CTSMI",0.,0);
	cpgswin(tmin,tmax,fmin,fmax);
	cpglab("Sub Integration","Channel","");
	cpgmtxt("R",3.0,0.5,0.5,"Frequency (MHz)");
	plot_image(ds,dynamic_levels);
	cpgbox("BCTSNI",0.,0,"BTSNI",0.,0);
      } else if (plot==wt_plot) {
	cpgeras();
	cpgswin(tmin,tmax,freqmin,freqmax);
	cpgbox("BCTSNI",0.,0,"CTSMI",0.,0);
	cpgswin(tmin,tmax,fmin,fmax);
	cpglab("Sub Integration","Channel","");
	cpgmtxt("R",3.0,0.5,0.5,"Frequency (MHz)");
	plot_image(wt,dynamic_levels);
	cpgbox("BCTSNI",0.,0,"BTSNI",0.,0);
      } else if (plot==prof_plot) {
	cpgeras();
	plot_profile(fp);
      }
    }

    // Mark zapped profiles
    for (unsigned i=izap;i<zap_list.size();i++) {
      cpgsci(2);
      if (plot==time_plot) {
	cpgrect(0.0,1.0,zap_list[i].sub0+0.5,zap_list[i].sub1+0.5);
      } else if (plot==freq_plot) {
	cpgrect(0.0,1.0,zap_list[i].chan0+0.5,zap_list[i].chan1+0.5);
      } else if (plot==both_plot || plot==wt_plot) {
	cpgrect(zap_list[i].sub0+0.5,zap_list[i].sub1+0.5,zap_list[i].chan0+0.5,zap_list[i].chan1+0.5);
      }
      cpgsci(1);
    }

    // Read cursor position
    if (click==false) {
      if (curs==horizontal) mode=5;
      if (curs==vertical) mode=6;
      if (curs==both) mode=7;
      cpgband(mode,0,0,0,&x0,&y0,&ch);
    } else if (click==true) {
      if (curs==horizontal) mode=3;
      if (curs==vertical) mode=4;
      if (curs==both) mode=2;
      cpgband(mode,0,x0,y0,&x1,&y1,&ch);
    }

    // Time mode
    if (ch==CMD_TIME) {
      plot=time_plot;
      click=true;
      redraw=true;
      replot=true;
      click=false;
      curs=horizontal;
    }

    // Freq mode
    if (ch==CMD_FREQ) {
      plot=freq_plot;
      click=true;
      redraw=true;
      replot=true;
      click=false;
      curs=horizontal;
    }

    // DS mode
    if (ch==CMD_BOTH) {
      plot=both_plot;
      click=true;
      redraw=true;
      replot=true;
      click=false;
    }

    // DS mode
    if (ch==CMD_WT) {
      plot=wt_plot;
      click=true;
      redraw=true;
      replot=true;
      click=false;
    }

    // Profile mode
    if (ch==CMD_PROF) {
      plot=prof_plot;
      click=true;
      redraw=true;
      replot=true;
      click=false;
    }
    
    // Left mouse click for zooming
    if (ch=='A') {
      if (click==false) {
	click=true;
	redraw=false;
	replot=false;
      } else if (click==true) {
	if (plot==time_plot) {
	  if (curs==horizontal) {
	    tmin=(y0<y1) ? y0 : y1;
	    tmax=(y0>y1) ? y0 : y1;
	  }
	} else if (plot==freq_plot) {
	  if (curs==horizontal) {
	    fmin=(y0<y1) ? y0 : y1;
	    fmax=(y0>y1) ? y0 : y1;
	  }
	} else if (plot==both_plot || plot==wt_plot) {
	  if (curs==horizontal || curs==both) {
	    fmin=(y0<y1) ? y0 : y1;
	    fmax=(y0>y1) ? y0 : y1;
	  }
	  if (curs==vertical || curs==both) {
	    tmin=(x0<x1) ? x0 : x1;
	    tmax=(x0>x1) ? x0 : x1;
	  }
	}
	redraw=true;
	click=false;
	replot=true;
      }
      continue;
    }

    // Unzoom
    if (ch==CMD_UNZOOM) {
      redraw=true;
      replot=true;
      unzoom=true;
      click=false;
      izap=zap_list.size();      
      continue;
    }

    // Levels
    if (ch==CMD_LEVELS) {
      if (dynamic_levels==true)
	dynamic_levels=false;
      else if (dynamic_levels==false)
	dynamic_levels=true;
      redraw=true;
      replot=true;
      click=false;
    }

    // Zap
    if (ch=='X') {
      if (click==false) { // Zap a single row or column
	if (plot==time_plot) {
	  zap.sub0=zap.sub1=(int) y0;
	  zap.type=time_zap;
	} else if (plot==freq_plot) {
	  zap.chan0=zap.chan1=(int) y0;
	  zap.type=freq_zap;
	} else if (plot==both_plot || plot==wt_plot) {
	  if (curs==both) {
	    zap.sub0=zap.sub1=(int) x0;
	    zap.chan0=zap.chan1=(int) y0;
	    zap.type=both_zap;
	  } else if (curs==horizontal) {
	    zap.chan0=zap.chan1=(int) y0;
	    zap.type=freq_zap;
	  } else if (curs==vertical) {
	    zap.sub0=zap.sub1=(int) x0;
	    zap.type=time_zap;
	  }
	}
      } else if (click==true) { // Zap a range
	if (plot==time_plot) {
	  zap.sub0=(int) y0;
	  zap.sub1=(int) y1;
	  zap.type=time_zap;
	} else if (plot==freq_plot) {
	  zap.chan0=(int) y0;
	  zap.chan1=(int) y1;
	  zap.type=freq_zap;
	} else if (plot==both_plot || plot==wt_plot) {
	  if (curs==both) {
	    zap.sub0=(int) x0;
	    zap.chan0=(int) y0;
	    zap.sub1=(int) x1;
	    zap.chan1=(int) y1;
	    zap.type=both_zap;
	  } else if (curs==horizontal) {
	    zap.chan0=(int) y0;
	    zap.chan1=(int) y1;
	    zap.type=freq_zap;
	  } else if (curs==vertical) {
	    zap.sub0=(int) x0;
	    zap.sub1=(int) x1;
	    zap.type=time_zap;
	  }
	}
      }
      apply_zap(arch,&zap);
      zap_list.push_back(zap);
      replot=true;
      redraw=true;
      click=false;
      continue;
    }
      
    // Middle mouse click = redraw diagnostics
    if (ch==CMD_UPDATE) {
      redraw=true;
      replot=true;
      izap=zap_list.size();
      click=false;
      continue;
    }

    // Toggle dedispersion
    if (ch==CMD_DEDISPERSE) {
      if (arch->get_dedispersed()) {
	// Reload clone
	arch = orig_arch->clone();

	// Reapply zap
	for (unsigned i=0;i<zap_list.size();i++) {
	  zap=zap_list[i];
	  apply_zap_range(arch,&zap);
	}
      } else {
	// Rededisperse
	arch->dedisperse();
      }
      // Remove baseline
      arch->remove_baseline();
      redraw=true;
      replot=true;
      click=false;
      continue;
    }

    // Undo last zap
    if (ch==CMD_UNDO) {
      if (!zap_list.empty()) {
	zap=zap_list.back();
	apply_zap_range(arch,&zap,true,orig_arch);
	zap_list.pop_back();
	// Reapply whole list in case of overlapping zaps
	for (unsigned i=0;i<zap_list.size();i++) {
	  zap=zap_list[i];
	  apply_zap_range(arch,&zap);
	}
      }
      redraw=true;
      replot=true;
      click=false;
      continue;
    }

    // Horizontal cursor
    if (ch==CMD_FREQMODE) {
      curs=horizontal;
      click=false;
      continue;
    }

    // Vertical cursor
    if (ch==CMD_TIMEMODE && plot==both) {
      curs=vertical;
      click=false;
      continue;
    }

    // both cursor
    if (ch==CMD_BOTHMODE && plot==both) {
      curs=both;
      click=false;
      continue;
    }

    // Generate PSRSH script
    if (ch==CMD_PSRSH) {
      FILE *file;
      file=fopen(psrsh_filename.c_str(),"w");
      fprintf(file,"#!/usr/bin/env psrsh\n\n# Run with psrsh -e <ext> <script.psh> <archive.ar>\n\n");

      // Full subint/channel zap
      for (unsigned i=0;i<zap_list.size();i++) {
	int chan0=zap_list[i].chan0;
	int chan1=zap_list[i].chan1;
	int sub0=zap_list[i].sub0;
	int sub1=zap_list[i].sub1;
	if (zap_list[i].type==freq_zap) {
	  if (chan0==chan1)
	    fprintf(file,"zap chan %d\n",chan0);
	  else
	    fprintf(file,"zap chan %d-%d\n",chan0,chan1);
	} else if (zap_list[i].type==time_zap) {
	  if (sub0==sub1) 
	    fprintf(file,"zap subint %d\n",sub0);
	  else
	    fprintf(file,"zap subint %d-%d\n",sub0,sub1);
	} else if (zap_list[i].type==both_zap) {
	  if (chan0==chan1 && sub0==sub1) {
	    fprintf(file,"zap such %d,%d\n",sub0,chan0);
	  } else {
	    fprintf(file,"zap such");
	    for (int sub=sub0;sub<=sub1;sub++) 
	      for (int chan=chan0;chan<=chan1;chan++) 
		fprintf(file," %d,%d",sub,chan);
	    fprintf(file,"\n");
	  }
	}
      }
      fclose(file);
      printf("PSRSH script written\n");
      click=false;
      continue;
    }
    
    // Save
    if (ch==CMD_SAVE) {
      // Apply zap to the original file
      for (unsigned i=0;i<zap_list.size();i++) {
	zap=zap_list[i];
	apply_zap_range(orig_arch,&zap);
      }
      orig_arch->unload(output_filename);
      click=false;
      printf("File saved\n");
      continue;
    }

    // Show help
    if (ch==CMD_HELP) {
      usage_interactive();
      click=false;
      continue;
    }
    
  } while (ch!=CMD_QUIT);


  cpgend();

  free(fp.z);
  free(tp.z);
  free(ds.z);

  return 0;
}

