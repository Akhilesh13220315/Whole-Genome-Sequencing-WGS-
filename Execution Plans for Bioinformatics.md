Project 1: Analysis of Protein Sequences to Predict Structure and Function
Objective
Analyze a protein sequence to predict its 3D structure (using AlphaFold) and function (e.g., domains, motifs via InterPro), enabling insights into biological roles.
Tools

Biopython (for sequence handling).
AlphaFold (via ColabFold or local install).
InterProScan (for function prediction).
PyMOL/VMD for visualization.

Reference Data

Use UniProt entry P00760 (Trypsin from Bos taurus) as an example sequence. Download FASTA: wget https://rest.uniprot.org/uniprotkb/P00760.fasta.
Open dataset: AlphaFold DB for predicted structures (e.g., AF-P00760-F1-model_v4.pdb from https://alphafold.ebi.ac.uk/entry/P00760).
PDB structure for validation: 1S0Q (Trypsin crystal structure) from RCSB PDB. Download: wget https://files.rcsb.org/download/1S0Q.pdb.

Steps

Sequence Retrieval: Use Biopython to fetch the sequence.

Code:
pythonfrom Bio import SeqIO
sequence = SeqIO.read("P00760.fasta", "fasta")
print(sequence.seq)



Multiple Sequence Alignment (MSA): Generate MSA for evolutionary context (input for AlphaFold).

Use MMseqs2 or ColabFold's built-in MSA.
Command: In ColabFold, upload FASTA and run.


Structure Prediction: Run AlphaFold to predict 3D structure.

Use Google Colab: Open https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb, input sequence, and predict.
Output: PDB file with predicted model.


Function Prediction: Annotate sequence with InterPro for domains/motifs.

Command: interproscan.sh -i P00760.fasta -f tsv (install via conda).
Parse TSV for functional annotations (e.g., serine protease domain).


Validation: Compare predicted structure to known PDB using RMSD (e.g., via PyMOL align command).
Analysis: Interpret results (e.g., active sites from InterPro).

Visualization

PyMOL: Load predicted PDB: pymol AF-P00760-F1-model_v4.pdb, color by secondary structure (util.rainbow), and show surface (show surface).
VMD: vmd 1S0Q.pdb, select protein (protein), draw cartoon representation (NewCartoon).
Tutorial: Follow PyMOL beginner guide for rendering images.

Expected Outputs

Predicted PDB file.
InterPro TSV with functional annotations.
RMSD value (e.g., <2Å for accurate prediction).
Visualized images (e.g., PNG exports from PyMOL).

Notes

Runtime: ~10-30 min on Colab GPU.
Limitations: AlphaFold excels for monomers; for complexes, use AlphaFold-Multimer.
Scaling: Batch process multiple sequences via AlphaFold DB API.

Project 2: Protein Docking
Objective
Predict binding poses and affinities between a protein and ligand (or another protein) to study interactions, e.g., for drug design.
Tools

AutoDock Vina/AutoDock Tools (ADT) for docking.
Biopython for file handling.
PyMOL for visualization.

Reference Data

Protein: PDB ID 1S0Q (Trypsin). Download as above.
Ligand: Benzamidine (inhibitor) from PubChem CID 2332 (SMILES: C1=CC=C(C=C1)C(=N)N). Download SDF: wget https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2332/record/sdf -O ligand.sdf.
Open dataset: Structural Protein Sequences from Kaggle (includes PDB files for docking benchmarks).

Steps

Prepare Protein and Ligand: Clean PDB (remove waters/heteroatoms) using ADT.

In ADT: Load 1S0Q.pdb, select protein, save as PDBQT.
Convert ligand SDF to PDBQT: mk_prepare_ligand.py -i ligand.sdf -o ligand.pdbqt (from AutoDock).


Define Grid Box: Use ADT to set search space around active site (e.g., center on benzamidine binding pocket).

Grid parameters: Size 20x20x20 Å, center at known site coordinates.


Run Docking: Perform docking with AutoDock Vina.

Command: vina --receptor protein.pdbqt --ligand ligand.pdbqt --center_x=10 --center_y=20 --center_z=30 --size_x=20 --size_y=20 --size_z=20 --out docked.pdbqt.


Analyze Results: Score poses by binding energy (lowest is best).

Use ADT to view top poses.


Validation: Compare to known complex (e.g., PDB 1S0Q has bound benzamidine).

Visualization

PyMOL: Load protein.pdb and docked.pdbqt, show sticks for ligand (show sticks, ligand), and hydrogen bonds (dist hbond, all, all).
VMD: Load files, highlight interactions with Graphics > Representations > Bonds.
Tutorial: AutoDock analysis in PyMOL.

Expected Outputs

Docked PDBQT files with poses.
Binding affinities (e.g., -7 kcal/mol).
Interaction maps (H-bonds, hydrophobic).

Notes

Runtime: ~5-15 min per docking.
For protein-protein: Use HADDOCK or ZDOCK instead.
Limitations: Rigid docking; for flexibility, use AutoDockFR.

Project 3: Homology Modeling
Objective
Build a 3D model of a target protein using a template with sequence similarity, useful when no experimental structure exists.
Tools

MODELLER for modeling.
BLAST for template search.
PyMOL for visualization.

Reference Data

Target sequence: UniProt P00760 (Trypsin, as above).
Template: PDB 1S0Q (high identity template).
Open dataset: InterPro for sequence annotations; download alignments from UniProt.

Steps

Template Search: BLAST target sequence against PDB.

Use NCBI BLAST: Input FASTA, select PDB database, choose top hit (e.g., 1S0Q).


Sequence Alignment: Align target and template.

Code (MODELLER Python):
pythonfrom modeller import *
env = Environ()
aln = Alignment(env)
mdl = Model(env, file='1S0Q', model_segment=('FIRST:A','LAST:A'))
aln.append_model(mdl, align_codes='1S0Q')
aln.append(file='P00760.ali', align_codes='target')
aln.align2d()
aln.write(file='alignment.ali', alignment_format='PIR')



Model Building: Generate models with MODELLER.

Script: Set a.starting_model=1; a.ending_model=5; a.make().


Model Assessment: Evaluate with DOPE score; select best model.
Refinement: Loop modeling if needed (loop.refine()).

Visualization

PyMOL: Load model.pdb, align to template (align model, 1S0Q), color by RMSD.
VMD: Compare structures with measure rmsd.
Tutorial: MODELLER with PyMOL superposition.

Expected Outputs

Modeled PDB files (e.g., target.B99990001.pdb).
Alignment PIR file.
DOPE scores for model quality.

Notes

Runtime: ~10 min for small proteins.
Limitations: Requires >30% sequence identity for accuracy.
Scaling: Use SWISS-MODEL web server for quick runs.

Project 4: Molecular Dynamics Simulation
Objective
Simulate protein motion over time to study stability, folding, or interactions.
Tools

GROMACS for simulation.
PyMOL/VMD for trajectory analysis.

Reference Data

Protein: PDB 1S0Q (Trypsin).
Open dataset: OpenProteinSet for MSAs if needed for setup. Or use GROMACS tutorial data (lysozyme in water).

Steps

System Preparation: Clean PDB, add hydrogens.

Command: pdb2gmx -f 1S0Q.pdb -o processed.gro -water spce (CHARMM force field).


Solvation and Ionization: Add water box and ions.

editconf -f processed.gro -o boxed.gro -bt cubic -d 1.0
solvate -cp boxed.gro -cs spc216.gro -o solv.gro
gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral


Energy Minimization: Relax system.

grompp -f minim.mdp -c ions.gro -p topol.top -o em.tpr
mdrun -v -deffnm em


Equilibration: NVT then NPT.

Similar grompp and mdrun for nvt.mdp and npt.mdp.


Production Run: Simulate dynamics.

grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
mdrun -deffnm md (e.g., 10 ns).


Analysis: RMSD, RMSF, etc.

gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns



Visualization

VMD: Load .gro and .xtc: vmd processed.gro md.xtc, animate trajectory, plot RMSD.
PyMOL: Import trajectory plugin for viewing.
Tutorial: GROMACS with VMD for dynamics visualization.

Expected Outputs

Trajectory files (.xtc).
Analysis plots (e.g., RMSD.xvg).
Energy logs.