### Regress

The .pf files in this directory are intended for testing out
various verisons of python using the script regression.py which 
is in py_progs.  

Note that in some cases, there may have been changes to the .pf
files that are used to create particular models.

An example of this, is that beginning with python_82a the treatment
of reflection is read in has changed.  

So for example, in 82a

Surface.reflection.or.absorption(0=no.rerad,1=high.albedo,2=thermalized.rerad)                    0

sets the reflection, whereas earlier (when only disk reflection is allowed it is set by

Disk.illumination.treatment(0=no.rerad,1=high.albedo,2=thermalized.rerad,3=analytic)   0

To handle this, there are likely variations in the .pf files here, from the ones that would
be produced from the same models if one had entered everything interactively.


Ideally, the history of these changes would be documented here:

82a - Since reflection in 82a refers both to the star and disk, the variable  

Surface.reflection is now always requested, where as in earlier versions of python it was 
not needed, in situations where there was no disk.  AS a result, this rdpar line has been added
to 1d_sn.pf and star.pf


