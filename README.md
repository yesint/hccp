# Hierarchical Clustering of Correlation Patterns (HCCP)

HCCP is the algorithm and software for dynamic domain identification in proteins. See the following papers for details:

1. [Yesylevskyy, S. O., V. N. Kharkyanen, and A. P. Demchenko. Hierarchical clustering of the correlation patterns: New method of domain identification in proteins. Biophys. Chem. 2006. 119:84-93.](https://www.sciencedirect.com/science/article/pii/S0301462205001523)
2. [Yesylevskyy, S. O., V. N. Kharkyanen, and A. P. Demchenko. Dynamic protein domains: identification, interdependence and stability. Biophys J, 2006, 91, 670-685.](https://www.sciencedirect.com/science/article/pii/S0006349506717660?via%3Dihub)

The code was written in 2006-2009 in Fortran90 and was never published being available by request only. In 2025 I resurracted it by fixing minor compatibility issues with modern `gfortran` and uploaded to GitHub.

Despite being very old and written in Fortran HCCP is rather user-friendly and easy to modify.

When using HCCP please always cite the papers above.

## Compiling
You need `gfortran` to compile. Just run `make`.

## Running hccp
```hccp < inp```
where inp is your input file

Use the `sample_input` file as a template. The comments there are detailed enough.
The input file is format-specific! Do not add or remove lines - it won't work.
Only standard PDB files can be read, no CHARMm pdbs and the like.

## Output
Standard output contains:
- Natural number of clusters (most probable for this protein)
- Second natural number of clusters (the second probable)
- Mean intra-domain correlation and
- Mean interdomain correlation for the natural number of clusters.
 NOTE: These features are not published yet.

 The most useful information is written to the file `color.tcl`.
 The human-readable lines begin with #. Typical output looks like:
```
 #-- Number of clusters  2 -------
 ...
 # (  1|   2: 108)
 # (  2| 109: 254)
 # (  1| 255: 284)
 # (  2| 285: 306)
      ^  ^    ^
      |  |    End of segment
      |  Start of segment
      Domain
```
 Domains can consist of many segments, which are continuos in sequence.
 Chain identifiers are added if needed.

 Such blocks of information are written for ALL hierarchical levels.
 It is a good idea to read the file from the end - largest clusters are there.

 File `color.tcl` is re-written each time you run hccp with another protein. 
 Be careful not to overwrite your data!

## Using VMD to visualize domains
- Start VMD
- Type `cd <path to the directory where your files reside>` in console
- Load your pdb
- Type `source color.tcl`
 
 Now you can use the following commands:

- `next`    - show next hierarchical level.
- `prev`    - show previous hierarchical level.
- `last`    - show last level (2 domains). 
- `goto <N>`  - go to the level number N.

 The clusters are color-coded. Segments, which appear in blue are not assigned
 to any cluster larger then 5 reidues (they still can be assigned to smaller
 clusters).


  ## Possible problems
  - _PDB can not be read_: Check if this PDB is valid. Try to remove the ligands and hetero atoms. 
     Check the histidine residues - they should be marked HIS, other names 
     are not supported. 

  - _Incorrect number of residues is detected_:
   The number of residues is determined as a number of CA atoms. If you have
     some CA missed, you'll got fewer residues. Non-standard residues are also
     ignored, even if they have CA.

  - _Other problems:_ Please open GitHub issue
