#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <complex.h>
#include <fftw3.h>
#include <sndfile.h>

#include <math.h>

#include "gnuplot_i.h"

#define	FRAME_SIZE 1024
#define HOP_SIZE 1024

#define OCTAVE 12

static gnuplot_ctrl *h;
static fftw_plan plan;

static void
print_usage (char *progname)
{	printf ("\nUsage : %s <input file> \n", progname) ;
	puts ("\n"
		) ;

} 
static void
fill_buffer(double *buffer, double *new_buffer)
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

static int
read_n_samples (SNDFILE * infile, double * buffer, int channels, int n)
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

static int
read_samples (SNDFILE * infile, double * buffer, int channels)
{
  return read_n_samples (infile, buffer, channels, HOP_SIZE);
}

void
fft_init (complex in[FRAME_SIZE], complex spec[FRAME_SIZE])
{
  plan = fftw_plan_dft_1d(FRAME_SIZE, in, spec, FFTW_FORWARD, FFTW_ESTIMATE);
}

void
fft_exit (void)
{
  fftw_destroy_plan(plan);
}

void
fft_process (void)
{
  fftw_execute(plan);
}

char*
get_note_octave(int pitch){
	int note = pitch%12;
	switch(note){
		case 0:
			return "C";
		case 1:
			return "C*";
		case 2:
			return "D";
		case 3:
			return "D*";
		case 4:
			return "E";
		case 5:
			return "F";
		case 6:
			return "F*";
		case 7:
			return "G";
		case 8:
			return "G*";
		case 9:
			return "A";
		case 10:
			return "A*";
		case 11:
			return "B";
		default:
			return "error";
	}
}

int
main (int argc, char * argv [])
{	char 		*progname, *infilename;
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

	double decalage[FRAME_SIZE];

		
	while (read_samples (infile, new_buffer, sfinfo.channels)==1)
	  {
	    /* Process Samples */
	    printf("Processing frame %d\n",nb_frames);

	    /* hop size */
	    fill_buffer(buffer, new_buffer);

		for(int tho=0; tho<FRAME_SIZE; tho++){
			decalage[tho] = 0;
			for(int n=0; n<FRAME_SIZE-tho; n++){
				decalage[tho] += buffer[n]*buffer[n+tho];
			}
			decalage[tho] /= FRAME_SIZE;
		}

		double max = 0;
		int indice_max;
		for(int i=10; i<FRAME_SIZE; i++){
			if(decalage[i] > max){
				max = decalage[i];
				indice_max = i;
			}
		}

		printf("indice_max: %d\n", indice_max);

		double periode = (double)indice_max/44100;
		double freq_max = 1.0/periode;

		printf("freq_max : %f\n", freq_max);

	    
	    /* note */
	    int note = (int) round(57 + OCTAVE*log2(freq_max/440));
	    printf("note %d \n", note);
		char* note_octave = get_note_octave(note);
		printf("note_octave %s \n", note_octave);
	    

	    
	    
	    /* plot amplitude */
	    gnuplot_resetplot(h);
	    gnuplot_plot_x(h,decalage,FRAME_SIZE,"decalage");
	    sleep(1);
    
	    /* PLOT */
	    //gnuplot_resetplot(h);
	    //gnuplot_plot_x(h,buffer,FRAME_SIZE,"temporal frame");
	    //sleep(1);
    
    	    
	    nb_frames++;
	  }

	sf_close (infile) ;

	/* FFT exit */
	fft_exit();
	
	return 0 ;
} /* main */

