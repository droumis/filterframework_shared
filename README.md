**-------- flab analysis repo README----------------**

*See repo Wiki for analysis-specific documentation*

# **Branch management:** #
master <-- develop <-- <feature branch>

**General Repo Rules**

 **DO NOT** put anything into the shared develop branch with initials appended

* The develop branch is for basic, lab-shareable code
* Create a user-specific directory for yourself if you wish to put your kk_sj_sk_Kodez_useme_old_v4.m files

 For substantial modifications to develop, create a temporary working branch (ex: linearization-dev). 


**Codebase map:**

 **DFAnalysis:** analysis-related functions (‘dfa_...’, other)

* Example: dfa_calcriptrigspectrogram

 **DFScripts:** analysis-related scripts (‘dfs_...’)

* Example: dfs_riptriggeredspectrogram

 **DFFunctions:** filter framework data processing-related functions (‘get...’, ‘set…’, Iterators)

* Examples: get2dstate, getconstimes, eegprocess, geteegtimes, getvalideegtimes

 **Process:** functions and scripts that process data into the filter framework format

* Examples: createrewardinfo, createcellinfo, lineardayprocess



**General Info:**

 Where applicable, data is stored as nested matlab cell arrays.
 * {animal}{day}{epoch}{tetrode}{cell}