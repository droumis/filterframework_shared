**-------- flab analysis repo README----------------**

*See repo Wiki for analysis-specific documentation*

# **Branch management:** #
master <-- develop <-- <feature branch>

**General Repo Rules**

1. **DO NOT** put anything into the shared develop branch with initials appended

* The develop branch is for basic, lab-shareable code
* Create a user-specific directory for yourself if you wish to put your kk_sj_sk_Kodez_useme_old_v4.m files

2. For substantial modifications to develop, create a temporary working branch (ex: linearization-dev). 


**Codebase map:**

1. **DFAnalysis:** analysis-related functions (‘dfa_...’, other)

* Example: dfa_calcriptrigspectrogram

2. **DFScripts:** analysis-related scripts (‘dfs_...’)

* Example: dfs_riptriggeredspectrogram

3. **DFFunctions:** filter framework data processing-related functions (‘get...’, ‘set…’, Iterators)

* Examples: get2dstate, getconstimes, eegprocess, geteegtimes, getvalideegtimes

4. **Process:** functions and scripts that process data into the filter framework format

* Examples: createrewardinfo, createcellinfo, lineardayprocess



**General Info:**

1. Where applicable, data is stored as nested matlab cell arrays.
* {animal}{day}{epoch}{tetrode}{cell}