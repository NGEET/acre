# Answer Changing Regression Evaluator (ACRE)

A python toolset that evaluates land-model output (a la ALM/CLM) for regression testing

The `acre_driver.py` script is intended to diagnose the output of several single site runs, as a rapid visual and textual analysis on ecosystem response over a period of 5 - 50 years.  As a default, the user must provide the file-path to restart files. As a default, plotting is turned off.

The `acre_gridcomp.py` script is intended to diagnose the output of gridded CLM/ELM/FATES simulations, for rapid visualization.

Options for both include:
1. the ability to do a regression against another set of files (base).
2. analysis of history files (which have a 1-hour suggested output) (NOT `acre_gridcomp.py`?)
3. plotting (most of the analysis right now is visual, not much
                             really happens right now without plotting)


## Next capabilities that are in queue
- add more restart output diagnostics
- come up with a list of PASS/FAIL diagnostics
- h1,h2,h3 files
- add more checks
- optimize math and averaging

## FYI: REALLY HELPFULL DEBUG TOOL
UNcomment the imported "code" library and add the following call in the code where you want to add a debug stop: `code.interact(local=locals())`

## Example usage

For usage: `python acre_driver.py -h`

### Single site

`python acre_driver.py --plotmode --test-hist-pref=<fullpath> --base-hist-pref=<fullpath> --test-name=<test> --base-name=<text> --eval-id=<id-string>`

