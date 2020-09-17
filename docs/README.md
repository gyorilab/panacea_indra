Automated assembly of cell-type-specific mechanisms for the regulation of pain and inflammation with INDRA
==========================================================================================================

<img align="left" src="https://raw.githubusercontent.com/sorgerlab/indra/master/doc/indra_logo.png" width="134" height="100" />
INDRA integrates multiple text-mining systems and pathway databases to
automatically extract mechanistic knowledge from the biomedical literature and
through a process of knowledge assembly, build
executable models and causal networks. Based on profiling and perturbational
data, these models can be contextualized to be cell-type specific and used to
explain experimental observations or to make predictions.



In the context of the [`DARPA Panacea program`](https://www.darpa.mil/program/panacea),
the [`INDRA team`](https://indralab.github.io/) at the
[`Laboratory of Systems Pharmacology, Harvard Medical School`](https://hits.harvard.edu/)
is working on understanding the regulation of pain and inflammation with the
goal of finding new therapeutics using INDRA.

- [Ion-channel mechanism knowledge base](#ion-channel-kb)
- [Ion-channel inhibitor search](#ion-channel-inhibitor-search)
- [A self-updating model of pain mechanisms](#pain-machine)
- [Neuro-immune interactome](#neuroimmune)
- [General applications](#general)
- [Funding](#funding)

Ion-channel mechanism knowledge base
------------------------------------

<img align="left" src="https://raw.githubusercontent.com/indralab/panacea_indra/website/docs/ion_channel_kb_network.png" width="150" height="112" />

We used INDRA to assemble all mechanisms that 65 ion-channels that are
particularly important for nociception are involved in.

Each ion-channel's interactions can be browsed as networks at
[`NDEx`](http://ndexbio.org/#/networkset/8f22b3bf-21d8-11ea-bb65-0ac135e8bacf).

And the literature evidence can be inspected on the pages linked to below. These
pages alsow support expert curation of statements which are fed back to INDRA
to improve the models (see curation tutorial
[`here`](https://indra.readthedocs.io/en/latest/tutorials/html_curation.html).

[`CACNA1A`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/CACNA1A.html)
[`CACNA1B`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/CACNA1B.html)
[`CACNA1C`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/CACNA1C.html)
[`CACNA1D`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/CACNA1D.html)
[`CACNA1E`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/CACNA1E.html)
[`CACNA1F`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/CACNA1F.html)
[`CACNA1G`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/CACNA1G.html)
[`CACNA1H`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/CACNA1H.html)
[`CACNA1I`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/CACNA1I.html)
[`CACNA1S`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/CACNA1S.html)
[`HCN1`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/HCN1.html)
[`HCN4`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/HCN4.html)
[`KCNA1`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNA1.html)
[`KCNA10`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNA10.html)
[`KCNA2`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNA2.html)
[`KCNA3`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNA3.html)
[`KCNA4`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNA4.html)
[`KCNA5`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNA5.html)
[`KCNA6`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNA6.html)
[`KCNA7`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNA7.html)
[`KCNB1`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNB1.html)
[`KCNB2`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNB2.html)
[`KCNC1`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNC1.html)
[`KCNC2`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNC2.html)
[`KCNC3`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNC3.html)
[`KCNC4`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNC4.html)
[`KCND1`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCND1.html)
[`KCND2`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCND2.html)
[`KCND3`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCND3.html)
[`KCNF1`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNF1.html)
[`KCNG1`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNG1.html)
[`KCNG2`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNG2.html)
[`KCNG3`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNG3.html)
[`KCNG4`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNG4.html)
[`KCNH1`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNH1.html)
[`KCNH2`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNH2.html)
[`KCNH3`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNH3.html)
[`KCNH4`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNH4.html)
[`KCNH5`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNH5.html)
[`KCNH6`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNH6.html)
[`KCNH7`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNH7.html)
[`KCNH8`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNH8.html)
[`KCNK18`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNK18.html)
[`KCNK2`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNK2.html)
[`KCNMA1`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNMA1.html)
[`KCNN1`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNN1.html)
[`KCNQ1`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNQ1.html)
[`KCNQ2`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNQ2.html)
[`KCNQ3`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNQ3.html)
[`KCNQ4`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNQ4.html)
[`KCNQ5`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNQ5.html)
[`KCNS1`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNS1.html)
[`KCNS2`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNS2.html)
[`KCNS3`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNS3.html)
[`KCNV1`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNV1.html)
[`KCNV2`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/KCNV2.html)
[`SCN10A`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/SCN10A.html)
[`SCN11A`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/SCN11A.html)
[`SCN1A`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/SCN1A.html)
[`SCN2A`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/SCN2A.html)
[`SCN3A`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/SCN3A.html)
[`SCN4A`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/SCN4A.html)
[`SCN5A`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/SCN5A.html)
[`SCN8A`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/SCN8A.html)
[`SCN9A`](https://bigmech.s3.amazonaws.com/panacea/ion_channels/SCN9A.html)

Ion-channel inhibitor search
----------------------------

<img align="left" src="https://raw.githubusercontent.com/indralab/panacea_indra/website/docs/ion_channel_inhibitor_search.png" width="150" height="112" />
The ion-channel inhibitor search is an application that allows selecting some
ion-channels as "desirable" drug targets and others as "undesirable" drug targets.
The application then searches all INDRA Statements in the ion-channel knowledge
base, combined with the [`Small Molecule Suite / Target Affintiy Spectrum`](https://labsyspharm.shinyapps.io/smallmoleculesuite/)
data.

The ion-channel inhibitor search is available [`here`](http://34.230.33.149:5000/).

A self-updating model of pain mechanisms
----------------------------------------

<img align="left" src="https://raw.githubusercontent.com/indralab/panacea_indra/website/docs/painmachine_image.png" width="150" height="112" />
EMMAA (Ecosystem of Machine-maintained Models with Automated Analysis) makes
available a set of
computational models that are kept up-to-date using automated machine reading,
knowledge-assembly, and model generation, integrating new discoveries
immediately as they become available.

The EMMAA model representing pain mechanisms can be found
[`here`](https://emmaa.indra.bio/dashboard/painmachine).

Neuro-immune interactome
------------------------
Using transcriptional profiling of immune cells and neurons following injury, we are building
cell-type specific models of how immune cells and pain sensing neurons interact.

The latest set of cell-type specific interactions can be browsed on these pages:
[`DCs`](https://bigmech.s3.amazonaws.com/panacea/neuroimmune/dcs.html)
[`Monocytes`](https://bigmech.s3.amazonaws.com/panacea/neuroimmune/monocytes.html)
[`Dermal macrophages`](https://bigmech.s3.amazonaws.com/panacea/neuroimmune/dermal_macs.html)
[`Resident macrophages`](https://bigmech.s3.amazonaws.com/panacea/neuroimmune/resident_macs.html)
[`M2a macrophages`](https://bigmech.s3.amazonaws.com/panacea/neuroimmune/m2a.html)
[`M2b macrophages`](https://bigmech.s3.amazonaws.com/panacea/neuroimmune/m2b.html)

General INDRA applications
--------------------------
There are several applications built on top of INDRA that that are generally
applicable to biomedical research and can therefore also be used to study pain
mechanisms.
- [`INDRA database`](https://db.indra.bio): The INDRA database website provides
  a search interface to find INDRA Statements assembled from the biomedical
  literature, browse their supporting evidence, and curate any errors.  Az
  example search relevant to Panacea is Object: SCN10A to find entities that
  regulate the Nav 1.8 ion channel.
- [`INDRA network search`](https://network.indra.bio): The INDRA network search
  allows finding causal paths, shared regulators, and common targets between
  two entities. An example search relevant to Panacea is Subject: bradykinin,
  Object: PKC (see [`here`](https://network.indra.bio/query?query=2094257329)).
- [`Dialogue.bio`](http://dialogue.bio): The dialogue.bio website allows
  launching dedicated human-machine dialogue sessions where you can upload your
  data (e.g., DE gene lists or gene expression profiles), discuss relevant
  mechanisms, and build model hypotheses using simple English dialogue.
  For instance, you could try the following series of questions:
  "what is Nav1.8?", "what does it regulate?",
  "which of those are transcription factors?".

Funding
-------
This project is funded under the DARPA Panacea program (HR00111920022).
