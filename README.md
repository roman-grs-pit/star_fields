Produce a simulated field of stars, with approximately correct density, k-band luminosity function, and spectral type distrbutions.

After a git clone, some data needs to be downloaded and environoment variables need to be defined:

export star_fields=<where you git cloned from>/star_fields
export star_data=$star_fields/data
export PYTHONPATH=$PYTHONPATH:$star_fields/py

Download extinction maps: https://www.dropbox.com/scl/fi/ivq5qi5cjhv7tahxiqdl8/NICERmaps.tar?rlkey=s0ymcsilvt3wvnpr8ehpdzd65&dl=0
unpack (e.g., tar -xvf NICERmaps.tar)
export allskyext=<where you downloaded NICERmaps.tar>/NICERmaps

Then, the basic code should work, e.g.:

python
>>>import insertstars_new as insert
>>>insert.insertstars(0,0) #creates a text file with distribution of stars centered at ra,dec=0,0 within area of one Roman field of view

