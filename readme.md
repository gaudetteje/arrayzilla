# Arrayzilla - Brown University Microphone Array

Arrayzilla is an ultrasonic microphone array panel with 224 channels.  It simultaneously captures and records the echolocation signals of bats during behavioral experiments while a) stationary on a platform, or b) in flight approaching the array.

This documentation gives an overview of the Arrayzilla data processing functions.  Note that `az_process_data()` automates all of these steps for a large set of data.  Below are the descriptions for each function that is called by `az_process_data()`.

## Processing and data analysis functions 

Arrayzilla generates two data files for each side of the array. This is because the recorder hardware is mirrored on each side.

### Data preconditioning
- `az_align_data(fname1)` - validates and fixes data misalignment in 1st file
- `az_align_data(fname2)` - validates and fixes data misalignment in 2nd file

### Map data channels to array coordinates
- `array = az_define_array(arrayfile)` - generates the array coordinates, which are hard coded and must be updated manually, if changed

### Detect call events
- `az_detect_events(fname1,fname2)` - detects discontinuous recording events using header field information.  Discontinuous events may be due to explicitly triggered data, power cycling the digital boards on the array, not erasing previously recorded data, or data corruption and/or hardware reset.
- `az_split_events(fname1,fname2,events)` - detects echolocation calls in each event by searching through time series data

### Process each detected call to generate the beam pattern and bulk parameters
- `az_process_beams(fname1,fname2,call,array)` - takes each detected call event and array definition and processes the raw data.  This is where the bulk of the work is done in reconstructing the beam pattern.

### Plotting and visualization functions
- `plotArrayPoints(array)` - plots geometry of array elements and localized source position.  Angular coordinates are mapped onto a sphere around the source position.
- `plotTimeSeries(fname1,fname2,events)` - plots the time series data from several representative channels.  events may be a single event or call structure or multiple events or call structures.
- `plotAuxChan(hdr,events)` - plots the time series data on the auxiliary channel for all events.  events may be a single event or call structure or multiple events or call structures.
- `plotBeamPattern(beam{idx})` - plots the processed beam pattern with a variety of options (see help documentation)

### Auxiliary processing functions (optional)
- `genAudioTrack(ref,wavfile)` - creates a WAV file using the data from the strongest reference channel and appropriate timing over the entire recording
- `genVideoTrack(beam,ref,aviname,frequency)` - creates an AVI video (no audio) with the beam pattern plotted at a specified frequency


## References

For additional documentation on the array hardware and acoustic signal processing, see:

[1] Gaudette, J. E., Kloepper, L. N., Warnecke, M., & Simmons, J. A. (2014). High resolution acoustic measurement system and beam pattern reconstruction method for bat echolocation emissions. The Journal of the Acoustical Society of America, 135(1), 513–520. http://doi.org/10.1121/1.4829661

[2] Gaudette, J. E., & Simmons, J. A. (2014, August 22). Observing the Invisible: Using Microphone Arrays to Study Bat Echolocation. Acoustics Today, 16–23.

[3] DiCecco, J., Gaudette, J. E., & Simmons, J. A. (2013). Multi-component separation and analysis of bat echolocation calls. The Journal of the Acoustical Society of America, 133(1), 538–546. http://doi.org/10.1121/1.4768877

_Last updated: 10 March 2013_
