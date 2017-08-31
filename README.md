**-------- flab analysis repo README----------------**

*See repo Wiki for a detailed overview of the filter framework structure and data formats, as well as analysis-specific documentation*

**Branch management:**
master <-- develop <-- <feature branch>

**Data storage:**
 Where applicable, data is stored as nested matlab cell arrays.

 * {animal}{day}{epoch}{tetrode}{cell}

#**General Repo Rules**#

 **DO NOT** put anything into the shared develop branch with initials appended

* The develop branch is for basic, lab-shareable code

 For substantial modifications to develop, create a temporary working branch (ex: linearization-dev). 


#**Codebase map:**#

 **DFAnalysis:** analysis-related functions (‘dfa_...’, other)

* Example: dfa_calcriptrigspectrogram

 **DFScripts:** analysis-related scripts (‘dfs_...’)

* Example: dfs_riptriggeredspectrogram

 **DFFunctions:** filter framework data processing-related functions (‘get...’, ‘set…’, Iterators)

* Examples: get2dstate, getconstimes, eegprocess, geteegtimes, getvalideegtimes

 **Process:** functions and scripts that process data into the filter framework format

* Examples: createrewardinfo, createcellinfo, lineardayprocess

**Scripts** and **Functions:**  Plotting scripts and functions that exist outside of filter framework

* Example: plotraster