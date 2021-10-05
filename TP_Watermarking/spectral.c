#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>

#include <sndfile.h>
#include <fftw3.h>
#include <time.h>

#include <math.h>
#include <complex.h>

#include "gnuplot_i.h"

/* taille de la fenetre */
#define	FRAME_SIZE 2205 //44100*0.05
/* avancement */
#define HOP_SIZE 2205

static int echelle = 44100/FRAME_SIZE; 

static gnuplot_ctrl *h;
static fftw_plan plan;

static void print_usage (char *progname)
{	printf ("\nUsage : %s <input file> \n", progname) ;
	puts ("\n"
		) ;

} 

static void fill_buffer(double *buffer, double *new_buffer)
{
  int i;
  double tmp[FRAME_SIZE-HOP_SIZE];
  
  /* save */
  for (i=0;i<FRAME_SIZE-HOP_SIZE;i++)
    tmp[i] = buffer[i+HOP_SIZE];
  
  /* save offset */
  for (i=0;i<(FRAME_SIZE-HOP_SIZE);i++)
    {
      buffer[i] = tmp[i];
    }
  
  for (i=0;i<HOP_SIZE;i++)
    {
      buffer[FRAME_SIZE-HOP_SIZE+i] = new_buffer[i];
    }
}

static int read_n_samples (SNDFILE * infile, double * buffer, int channels, int n)
{

  if (channels == 1)
    {
      /* MONO */
      int readcount ;

      readcount = sf_readf_double (infile, buffer, n);

      return readcount==n;
    }
  else if (channels == 2)
    {
      /* STEREO */
      double buf [2 * n] ;
      int readcount, k ;
      readcount = sf_readf_double (infile, buf, n);
      for (k = 0 ; k < readcount ; k++)
	buffer[k] = (buf [k * 2]+buf [k * 2+1])/2.0 ;

      return readcount==n;
    }
  else
    {
      /* FORMAT ERROR */
      printf ("Channel format error.\n");
    }
  
  return 0;
} 

static int read_samples (SNDFILE * infile, double * buffer, int channels)
{
  return read_n_samples (infile, buffer, channels, HOP_SIZE);
}



//----------FFT----------
void fft_init(double complex *data_in, double complex *data_out){
	plan = fftw_plan_dft_1d(FRAME_SIZE, data_in, data_out, FFTW_FORWARD, FFTW_ESTIMATE);
}

void fft_exit(fftw_plan plan){
	fftw_destroy_plan(plan);
}


double get_indice_max(double *amplitude){
	double max = 0;
	double curr;
	int indice;

	for(int i=18000/echelle; i<22000/echelle; i++){
		curr = amplitude[i];
		if(curr > amplitude[i-1] && curr > amplitude[i+1]){
			if(curr > max){
				max = amplitude[i];
				indice = i;
			}
		}
	}
	return indice;
}

double interpolation_parabolique(int echantillon, double *amplitude){
	double amp_l = 20*log10(amplitude[echantillon-1]);
	double amp_c = 20*log10(amplitude[echantillon]);
	double amp_r = 20*log10(amplitude[echantillon+1]);

	double d = 0.5 * (amp_l - amp_r)/(amp_l - 2*amp_c + amp_r);
	return echantillon + d;
}

bool hidden_freq(double *amplitude){
	for(int i=18000/echelle; i<22000/echelle; i++){
		if(amplitude[i] > 10){
			return true;
		}
	}
	return false;
}

char get_event(double freq){
	int deltaA = fabs(freq - 19126);
	int deltaB = fabs(freq - 19585);
	int deltaC = fabs(freq - 20032);

	if(deltaA < deltaB	&&	deltaA < deltaC){
		return 'A';
	} else if(deltaB < deltaA	&&	deltaB < deltaC){
		return 'B';
	} else if(deltaC < deltaA	&&	deltaC < deltaB){
		return 'C';
	} else{
		return 'e';
	}
}

int main (int argc, char * argv []){
	char 		*progname, *infilename;
	SNDFILE	 	*infile = NULL ;
	SF_INFO	 	sfinfo ;

	progname = strrchr (argv [0], '/') ;
	progname = progname ? progname + 1 : argv [0] ;

	if (argc != 2)
	{	print_usage (progname) ;
		return 1 ;
		} ;

	infilename = argv [1] ;

	if ((infile = sf_open (infilename, SFM_READ, &sfinfo)) == NULL)
	{	printf ("Not able to open input file %s.\n", infilename) ;
		puts (sf_strerror (NULL)) ;
		return 1 ;
		} ;

	/* Read WAV */
	int nb_frames = 0;
	double new_buffer[HOP_SIZE];
	double buffer[FRAME_SIZE];

	/* Plot Init */
	h=gnuplot_init();
	gnuplot_setstyle(h, "lines");
	
	int i;
	for (i=0;i<(FRAME_SIZE/HOP_SIZE-1);i++)
	  {
	    if (read_samples (infile, new_buffer, sfinfo.channels)==1)
	      fill_buffer(buffer, new_buffer);
	    else
	      {
		printf("not enough samples !!\n");
		return 1;
	      }
	  }

	/* Info file */
	printf("sample rate %d\n", sfinfo.samplerate);
	printf("channels %d\n", sfinfo.channels);
	printf("size %d\n", (int)sfinfo.frames);

	double *amplitude = (double*) calloc(FRAME_SIZE, sizeof(double));
	double complex *data_in = (double complex*) calloc(FRAME_SIZE, sizeof(double complex));
    double complex *data_out = (double complex*) calloc(FRAME_SIZE, sizeof(double complex));

	fft_init(data_in, data_out);

	int frameA = 8*echelle;
	int frameB = 16*echelle;
	int frameC = 24*echelle;
	
	int max_indice;
	char event;
	char eventList[20];
	double freq[3];
	int j=0;
	float time[20];
	float event_time;
	int nb_event = 0;
	double h_freq, indice;

	while (read_samples (infile, new_buffer, sfinfo.channels)==1)
	{
	    /* Process Samples */
	    //printf("Processing frame %d\n",nb_frames);

	    /* hop size */
	    fill_buffer(buffer, new_buffer);

		for(int i=0; i<FRAME_SIZE; i++){
			data_in[i] = buffer[i]*(0.5-0.5*cos(2.0*M_PI*i/FRAME_SIZE)); //Hann
		}

		fftw_execute(plan);

		for(int i=0; i<FRAME_SIZE; i++){
			amplitude[i] = cabs(data_out[i]);
		}

		//CALIBRAGE
	    if(nb_frames == frameA || nb_frames == frameB || nb_frames == frameC){
			max_indice = get_indice_max(amplitude);
			freq[j] = max_indice * echelle;
			j++;
		}


		if(hidden_freq(amplitude)){
			indice = interpolation_parabolique(get_indice_max(amplitude), amplitude);
			h_freq = indice * echelle;
			event_time = (float)nb_frames/echelle;
			if(event_time - time[nb_event-1] > 0.1){
				time[nb_event] = event_time;
				event = get_event(h_freq);
				eventList[nb_event] = event;
			}
			nb_event++;
		}

		//PLOT
		// gnuplot_resetplot(h);
		// gnuplot_plot_x(h, amplitude, FRAME_SIZE*0.5, "spectral");
		// sleep(1);
    
	    nb_frames++;
	}

	// for(int i=0; i<3; i++){
	// 	printf("freq[%d] : %f\n", i, freq[i]);
	// }

	for(int i=0; i<(sizeof(time)/sizeof(time[0])); i++){
		if(time[i] != 0){
			printf("Event %c at %fs\n", eventList[i], time[i]);
		}
	}
	fft_exit(plan);

	free(data_in);
	free(data_out);
	free(amplitude);

	sf_close (infile) ;

	return 0 ;
} /* main */

