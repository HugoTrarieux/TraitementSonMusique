# TraitementSonMusique
Some TP from Sound and Musique Processing class, in M1 at University of Bordeaux




### TP_Waterarking ###

The goal of this lab is to write functions to detect certain events from a recorded signal.
The events we are looking for are sinusoids with high frequencies so that they are difficult to hear for the users.
We assume that frequencies above 18000Hz are not audible and must be We assume that frequencies above 18000Hz
are not audible and should be considered as event triggering frequencies.

Write an analysis function that takes a wav file as a parameter, detects the events present,
and then displays the times and type of events detected. Apply the previous function to the other music
files contained in the Watermarking directory and indicate the sequence of events obtained.



### TP_TEL ###

Frequency dialing telephone keypads generate a signal composed of a sum of two sinusoids of constant amplitude,
whose frequencies vary according to the key pressed. This is the DTMF norm (Dual Tone Multi-Frequency).
The goal of this lab is to write functions to detect phone numbers from a recorded signal.

