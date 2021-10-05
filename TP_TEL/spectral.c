#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <sndfile.h>
#include <fftw3.h>
#include <time.h>

#include <math.h>
#include <complex.h>

#include "gnuplot_i.h"

/* taille de la fenetre */
#define	FRAME_SIZE 4410 //44100*0.1
/* avancement */
#define HOP_SIZE 4410

static gnuplot_ctrl *h;
static fftw_plan plan;
char* phoneNumber[12];

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


double* get_indice_max(double *amplitude){
	double max[2] = {0, 0};
	double curr = 0;

	for(int i=1; i<FRAME_SIZE*0.5; i++){
		curr = amplitude[i];

		if(curr > amplitude[i-1] && curr > amplitude[i+1]){

			if(max[0] <= max[1]){
				if(curr > max[0]){
					max[0] = i;
				}
			} else{
				if(curr > max[1]){
					max[1] = i;			
				}
			}

		}
	}
	return max;
}

double interpolation_parabolique(int echantillon, double *amplitude){
	double amp_l = 20*log10(amplitude[echantillon-1]);
	double amp_c = 20*log10(amplitude[echantillon]);
	double amp_r = 20*log10(amplitude[echantillon+1]);

	double d = 0.5 * (amp_l - amp_r)/(amp_l - 2*amp_c + amp_r);
	return echantillon + d;
}

char convert_freq_to_num(double freq1, double freq2){
	int f1 = round(freq1);
	int f2 = round(freq2);

	if(f1 == 697){
		if(f2 == 1209)
			return '1';
		if(f2 == 1336)
			return '2';
		if(f2 == 1477)
			return '3';
	} else if(f1 == 770){
		if(f2 == 1209)
			return '4';
		if(f2 == 1336)
			return '5';
		if(f2 == 1477)
			return '6';
	} else if(f1 == 852){
		if(f2 == 1209)
			return '7';
		if(f2 == 1336)
			return '8';
		if(f2 == 1477)
			return '9';
	} else if(f1 == 941){
		if(f2 == 1209)
			return '*';
		if(f2 == 1336)
			return '0';
		if(f2 == 1477)
			return '#';
	} else{
		return '.';
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
	double *phase = (double*) calloc(FRAME_SIZE, sizeof(double));
	double complex *spectrum = (double complex*) calloc(FRAME_SIZE, sizeof(double complex));

	double complex *data_in = (double complex*) calloc(FRAME_SIZE, sizeof(double complex));
    double complex *data_out = (double complex*) calloc(FRAME_SIZE, sizeof(double complex));

	fft_init(data_in, data_out);

	double* max_indice;
	double newIndice1, newIndice2, max_freq[2];
	char number, previousNum = 0;

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
		
		max_indice = get_indice_max(amplitude);
		newIndice1 = interpolation_parabolique(max_indice[0], amplitude);
		newIndice2 = interpolation_parabolique(max_indice[1], amplitude);
		max_freq[0] = newIndice1 * sfinfo.samplerate/FRAME_SIZE;
		max_freq[1] = newIndice2 * sfinfo.samplerate/FRAME_SIZE;

		number = convert_freq_to_num(max_freq[0], max_freq[1]);
		if(number != previousNum){
			previousNum = number;
			if(number != '.'){
				strncat(phoneNumber, &number, 1);
			}
		}

	    /* PLOT */
	    gnuplot_resetplot(h);
	    gnuplot_plot_x(h, amplitude, 500, "spectral");
	    sleep(1);
    
	    nb_frames++;
	}
	printf("phoneNumber: %s\n", phoneNumber);

	fft_exit(plan);

	free(data_in);
	free(data_out);
	free(amplitude);
	free(phase);
	free(spectrum);


	sf_close (infile) ;

	return 0 ;
} /* main */

