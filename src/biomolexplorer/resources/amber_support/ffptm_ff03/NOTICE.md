# Third-party Forcefield_PTM data

The 32 individual PREP files are extracted without changes to atom types, charges or internal coordinates from ffptm.in. Each corresponding frcmod file is unchanged from ct400556v_si_002.zip. Extraction adds only the standalone PREP header and STOP terminator.

Khoury et al., Forcefield_PTM: Ab Initio Charge and AMBER Forcefield Parameters for Frequently Occurring Post-Translational Modifications, JCTC (2013), DOI 10.1021/ct400556v. https://acs.figshare.com/articles/journal_contribution/2344720

Publisher license: CC BY-NC 4.0, https://creativecommons.org/licenses/by-nc/4.0/ . The noncommercial condition applies to these parameter data; application code retains its own license. Cite the original work. No endorsement is implied.

Profiles use protein ff03 throughout. Published residue charge states are recorded separately in manifest.json. Terminal variants, crosslinks, missing heavy-atom graphs and alternative protonation states require dedicated validated models; no chemistry is silently substituted.

Source naming errata: duplicate hydrogen labels in DA1 and SEN are made unique in extracted PREP files. The exact atom indices and name substitutions are recorded in manifest.json adaptations. No charge, atom type, coordinate or bond parameter is changed.
