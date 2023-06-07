This folder contains python scripts that make use of the astropy package *astroplan* to determine which early science targets are observable for a given day/time.

*early_science_planning.py* makes use of defualt astroplan constraints on airmass (>2), twilight (civil), moon separation (>45deg), and moon illumination (dark/grey/bright time based on halpha surface brightness)

*early_science_planning_customconstraint.py* replaces the moon separation and moon illumination constraints with a custom constraint which combines these two properties to evaluate the additional sky brightness due to the moon. The allowed brightness is still ranked based on target halpha surface brightness.

When running the scripts, the start and end dates must be specified as well as whether plots should be made for each target/night showing the constraint conditions for each 1 hour block throughout the night. An example of one of these plots is shown in "30Dor_constraint_plot.png". To run for the month of July with plotting off in ipython enter:

`%run early_science_planning.py 2023-07-01 2023-08-01 False`

Both scripts generate and excel file ("EarlyScienceObservability_standard/custom_constraint.xlsx") which contains a sheet for each date evaluated. These sheets have a column for every 1 hour block of the night and a list of all targets observable in that block. Observable here means that all constraints are met for that date/time.
