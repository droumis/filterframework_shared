**#-------- flab analysis repo README----------------**
# **Branch management:** #
master <-- develop <-- '<user>-dev'

* new members please create a new branch <user>-dev (e.g. 'demetris-dev') from develop and merge working functions/scripts from their branch to develop
* master branch will be reserved for very polished, documented code that would be shareable with non-franklab members.

#**Codebase map:**# 

1. **DFAnalysis:** analysis-related functions (‘dfa_...’, other)

* Example: dfa_calcriptrigspectrogram

2. **DFScripts:** analysis-related scripts (‘dfs_...’)

* Example: dfs_riptriggeredspectrogram

3. **DFFunctions:** filter framework processing-related functions (‘get...’, ‘set…’, other)

* Examples: get2dstate, getconstimes, eegprocess, geteegtimes, getvalideegtimes