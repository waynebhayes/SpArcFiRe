#!/usr/bin/env python
# coding: utf-8

# # Example notebook: Search data with CasJobs
# 
# SciServer Compute can talk to other components of SciServer through a series of <em>modules</em>, one for each component. This example notebook shows how to use the <strong><code>SciServer.CasJobs</code></strong> module to search through SciServer datasets from within your Python scripts. CasJobs is SciServer's data search tool, allowing you to select data from any of our big data datasets and/or your own uploaded datasets.
# 
# You are welcome (encouraged!) to copy these examples into another folder and modify them to meet your needs. You can use them as a starting point to create your own scripts. Please do not edit this notebook directly, because your edits may be overwritten if changes to the SciServer modules require changes to these example notebooks.
# 
# To run the example Python scripts in this notebook, click in any of the Code cells below (the ones with the gray backgrounds). Click the play button at the top of the window (just below the menubar) to run the script, or press Shift-Enter. When you run a cell, its output of will appear directly below the cell.

# ## Import modules
# 
# Like any Python modules, the SciServer modules must be imported before being used. The next code block first imports the SciServer modules you will need for this example notebook, then imports some other required modules. Comments in the code block explain what each module does. To learn how to import other modules, see the Python 3.5 import documentation (https://docs.python.org/3.5/reference/import.html), or the documentation of the module(s) you are trying to import.

# In[ ]:


import SciServer
from SciServer import CasJobs     # Communicate between SciServer Compute and CasJobs
print('Imported SciServer modules')

import pandas                                # data analysis tools
import numpy as np                           # numerical tools
from datetime import datetime, timedelta     # date and timestamp tools
from pprint import pprint                    # print human-readable output
print('Imported other needed modules')


# ## Get help
# 
# At any point after the modules are imported, you can type "help (<em>name of module</em>)" to read the documentation for that module. This is true for all SciServer modules and most other modules as well. Try it below.

# In[ ]:


# Read the help document for the entire SciServer package
# help(SciServer)


# In[ ]:


# Read the help document for the CasJobs module
# help(CasJobs)


# ### Optional but very helpful: create convenience functions
# 
# You can define your own functions in Python with the <code>def</code> (define) command. Once you define a function, you can use it later in the same notebook by calling its name and argument(s).
# 
# The Code cell below creates the following two "convenience functions" that make it easier to work with CasJobs output. These functions are used by Code cells later in this notebook. Run the Code cell to load these functions into memory; the cell will print a confirmation message. You can then use the two functions in later Code cells. Further examples in this notebook demonstrate how to use these functions.

# In[ ]:


# PYTHON CONVENIENCE FUNCTIONS USEFUL FOR WORKING WITH CASJOBS

def tables_formatted(tables):   # better formatted printing of a tables dictionary (output of get_tables)
# Returns the following information about the tables in your MyDB (as a Python dictionary object):
### Size: size of the table (in kB)
### Name: the name of the table
### Rows: the number of rows the table contains
### Date: the date of the table's creation, as the number of 10-microsecond intervals elapsed 1 AD

    import pandas
    from datetime import datetime
    
    tables = sorted(tables, key=lambda k: k['Name']) # alphabetize by table name
    
    for thistable in tables:
        print('Table name:\t',thistable['Name'])
        print('Rows:\t\t {:,.0f}'.format(thistable['Rows']))
        print('Size (kB):\t {:,.0f} '.format(thistable['Size']))

        cjCreateDate = thistable['Date']
        createsec = cjCreateDate / 10000000  # Divide by 10 million to get seconds elapsed since 1 AD
        firstday = datetime(1, 1, 1, 0, 0)   # Save 1 AD as "firstday"
        created = firstday + timedelta(seconds=createsec)  # Get calendar date on which table was created     
        print('Created time:\t',created.strftime('%Y-%m-%d %H:%M:%S'))
        print('\n')
        

def jobDescriber(jobDescription):
    # Prints the results of the CasJobs job status functions in a human-readable manner
    # Input: the python dictionary returned by getJobStatus(jobId) or waitForJob(jobId)
    # Output: prints the dictionary to screen with readable formatting
    import pandas
    
    if (jobDescription["Status"] == 0):
        status_word = 'Ready'
    elif (jobDescription["Status"] == 1):
        status_word = 'Started'
    elif (jobDescription["Status"] == 2):
        status_word = 'Cancelling'
    elif (jobDescription["Status"] == 3):
        status_word = 'Cancelled'
    elif (jobDescription["Status"] == 4):
        status_word = 'Failed'
    elif (jobDescription["Status"] == 5):
        status_word = 'Finished'
    else:
        status_word = 'Status not found!!!!!!!!!'

    print('JobID: ', jobDescription['JobID'])
    print('Status: ', status_word, ' (', jobDescription["Status"],')')
    print('Target (context being searched): ', jobDescription['Target'])
    print('Message: ', jobDescription['Message'])
    print('Created_Table: ', jobDescription['Created_Table'])
    print('Rows: ', jobDescription['Rows'])
    wait = pandas.to_datetime(jobDescription['TimeStart']) - pandas.to_datetime(jobDescription['TimeSubmit'])
    duration = pandas.to_datetime(jobDescription['TimeEnd']) - pandas.to_datetime(jobDescription['TimeStart'])
    print('Wait time: ',wait.seconds,' seconds')
    print('Query duration: ',duration.seconds, 'seconds')
        
print('Created functions')


# ## What data can I search?
# 
# CasJobs allows you to search many different datasets, referred to as <strong>contexts</strong> (they are known as contexts rather than databases, so that they can be described independently of the databases in which they are stored). Each context consists of one or more tables containing data or metadata related to a single aspect of the full dataset.

# ### Get a list of contexts
# 
# At the moment, the SciServer.CasJobs module does not have a function to list available contexts. The best way to see what contexts are available to you is to log in to <a href="http://skyserver.sdss.org/casjobs/" target="_blank">CasJobs</a> (link opens in a new window). Once you are logged in, you should see the Query page. Look for the <strong>Contexts</strong> dropdown menu toward the top left of the page, just above the big textbox. The values in that dropdown list show the contexts you can search, both directly in CasJobs and in Compute.

# ### Show data tables in a context
# 
# Once you know what context you want to search, you can use the <strong>CasJobs.getTables(context)</strong> function to show the data tables in that context. The Code cell below gives commands to list all tables in a context. Set the value of <em>this_context</em> to be the context you want to see. The function CasJobs.getTables(context) returns a list of Python dictionaries, one dictionary per table.
# 
# Each dictionary in the list contains the following information about one table:
# <ul>
# <li><em>Date:</em> the number of 10-millisecond intervals since the table was created</li>
# <li><em>Name:</em> the name of the table</li>
# <li><em>Rows:</em> the number of rows in the table</li>
# <li><em>Size:</em> the size of the table in kilobytes</li>
# </ul>
# 
# The code cell gives two options for printing the list of tables: using Python's <code>pprint</code> library or using the <code>tables_formatted(tableList)</code> convenience function defined above. The convenience function sorts the list of tables alphabetically by name, and displays the dates into datetime values. Try uncommenting and commenting those lines in the Code cell below to see what both options do.

# In[ ]:

#check_token = SciServer.Authentication.getToken()

my_token = SciServer.Authentication.login('portmanm','uB53!BYLQXr5XhN')
#this_context = "MyDB"    # Your MyDB
this_context = 'DR7'   # SDSS Data Release 14

tables = CasJobs.getTables(context=this_context)
#print('Tables in '+this_context+':\n')


#pprint(tables)   # Standard human-readable printing using Python's pprint module
# tables_formatted(tables)  # Sorting and better printing using a convenience function


# ## Run a quick query and get results directly
# 
# Now that you know what contexts (datasets) are available to you, and you know what tables can be found in those contexts, you are ready to write and submit a query to that context. A query is a request for data, written in SQL (Structured Query Language), a programming language designed for efficient database searches. 
# 
# SkyServer features a <a href="http://skyserver.sdss.org/public/en/help/howto/search/searchhowtohome.aspx" target="_blank">tutorial for learning SQL</a>, as well as <a href="http://skyserver.sdss.org/public/en/help/docs/sql_help.aspx" target="_blank">tips for writing good queries</a> and a long list of <a href="http://skyserver.sdss.org/public/en/help/docs/realquery.aspx" target="_blank">sample queries</a> that you can adapt to create your queries (links open in new windows).
# 
# Once you have written a query, you can get results by running (executing) it in CasJobs. To run a query in CasJobs directly from a Code cell in SciServer Compute, use the <strong><code>CasJobs.executeQuery(sql,...)</code></strong> function. The function takes as input a string containing a properly-formatted SQL query (and optional parameters listed below), and returns a table containing query results (in one of several formats with a default of a <a href="https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html" target="_blank">pandas dataframe</a>).
# 
# The <em>sql</em> parameter is required. The <em>context</em> parameter is recommended to explicitly state the context to which the query will be submitted; the default value is 'MyDB'. For a full list of parameters taken by the <code>CasJobs.executeQuery(sql,...)</code> function, see <a href="http://www.sciserver.org/docs/sciscript-python/SciServer.html#module-SciServer.CasJobs" target="_blank">its documentation on the SciServer documentation site</a>.
# 
# ### Example
# 
# The Code cell below gives an example of a query that will run quickly in CasJobs and return 10 results. It also includes (commented out) examples of an invalid SQL statement and a valid statement that returns no results, so you can see what happens in those cases.

# In[ ]:


# Example of the function CasJobs.executeQuery(sql, ...)
#   This function executes a quick SQL query and store results in a pandas dataframe.

# A valid query that returns results.
### Return IDs, positions, and g magnitudes (brightness measurements in red wavelength of light)...
#### for 10 galaxies with reliable imaging data.

myquery = 'select top 10 objid, ra, dec, r ' # note the space at the end of this string - important
myquery += 'from galaxy '
myquery += 'where clean = 1'

# A badly-formed query
#badquery = "ceci n''est pas une query"   # note substitution of '' for single quote mark

# A valid query that returns no results
#zeroquery = "select top 10 * from specobj where 0=1"

# Execute a query to the DR14 context, return a pandas dataframe, then show it
#df = CasJobs.executeQuery(sql=myquery, context=this_context)
# df


# ## Run a longer query and store results in your MyDB
# 
# The example above shows a quick query. Quick queries are limited to 60 seconds of processing time. For longer queries, CasJobs has a system of <strong>jobs</strong>. When you submit a query to CasJobs, the system creates a job for your query, with a unique <code>jobId</code>. CasJobs then runs this job in the background as server resources permit. When a job completes, it writes the query results into your MyDB personal database space, into a table that you specify. You can later retrieve and use the query results associated with the job by querying that table in your MyDB (that is, by executing a query with <code>context='MyDB'</code>).
# 
# ### Submitting a query as a job to CasJobs
# 
# To run a long query with MyDB output from SciServer Compute, you would use a different approach from the one you used in the previous Code cell. Submit the query using the function <strong><code>CasJobs.submitJob(sql,context)</code></strong>. The context used in this function is the context whose data you want to search (since the results will automatically be written into the <code>MyDB</code> context).
# 
# The <em>sql</em> parameter is required. The <em>context</em> parameter is highly recommended to explicitly state the context to which the query will be submitted; the default value is 'MyDB'. For a full list of parameters taken by the <code>CasJobs.submitJob(sql, context)</code> function, see <a href="http://www.sciserver.org/docs/sciscript-python/SciServer.html#module-SciServer.CasJobs" target="_blank">its documentation on the SciServer documentation site</a>.
# 
# Importantly, the function <strong><code>CasJobs.submitJob(sql,context)</code></strong> does NOT return query results as output. It returns only an integer <code>jobID</code> of the job that corresponds to your query. To retrieve the data and use it in a Compute script, you must still use get them from your MyDB by running <strong><code>CasJobs.executeQuery(sql,context='MyDB',...)</code></strong>, as described in the previous section.
# 
# 
# ### Check on running jobs
# 
# When you submit a job to CasJobs, it runs in the background. The time required to finish a job is hard to predict, because it can vary widely based on the efficiency of your query and the current load on the servers that power CasJobs &mdash; and the function <code>CasJobs.SubmitQuery(sql,context)</code> does not does not indicate when the job has completed, and thus when your results are ready in your MyDB.
# 
# We have created two functions to check on the status of running jobs and signal when those jobs are completed:
# 
# <ul>
# <li><strong><code>SciServer.CasJobs.waitForJob(jobId,...)</code></strong> takes as input a jobID (integer) and checks the status of that job every few seconds. It will display a "Waiting..." message while the job is running, and return its output once the job completes: a Python dictionary of metadata describing the job. Adding the optional argument <em>verbose=False</em> will suppress the "Waiting..." messages.</li>
# <li><strong><code>SciServer.CasJobs.waitForJob(jobId)</code></strong> also takes the jobID as input and returns the same type of metadata dictionary as output, but does not continuously check job status. Instead, it checks once and returns the status of that job at that moment. If the function returns output saying that the job has not finished, you will need to check again later by calling the function again.</li>
# </ul>
# 
# To view the metadata dictionary that results from of either of those functions, you can either use the imported <code>pprint</code> module, or use the <code>jobDescriber(result_dict)</code> convenience function defined above.
# 
# 
# ### Example
# 
# The Code cell below gives an example of a query that will likely be too long for the CasJobs quick queue. The script prints out a waiting message while the query runs, and then writes the query results into a table in your MyDB personal database space. The name of that table to be written is specified by the string variable <code>bigtablename</code> in the first line of the script. Be sure that this table does not already exist, or the query will fail with the message: <code>There is already an object named ... in the database.</code>
# 
# The last two lines of the Code cell below give you options for viewing a summary of the results of the job. You can use the <code>pprint</code> function from Python, or you can use the <code>jobDescriber(jobDescription)</code> convenience function defined above.

# In[ ]:


bigtablename = 'hugetable'

# Example of a longer query: get magnitudes and sizes (Petrosian radii) of one million galaxies
verylongquery = 'select top 10 objid, ra, dec \n'
verylongquery += 'u, g, r, i, z, err_u, err_g, err_r, err_i, err_z, petror90_r \n'
verylongquery += 'into mydb.' + bigtablename + '\n'
verylongquery += 'from galaxy\n'
verylongquery += 'where clean = 1'

print('Submitting query:\n',verylongquery)
print('\n')

thisjobid = CasJobs.submitJob(sql=verylongquery, context=this_context)

print('Job submitted with jobId = ',thisjobid)
print('\n')

waited = CasJobs.waitForJob(jobId=thisjobid)      # waited is a dummy variable; just print wait msg
jobDescription = CasJobs.getJobStatus(thisjobid)

print('\n')
print('Information about the job:')

#pprint(jobDescription)
jobDescriber(jobDescription)


# ## Thank you!
# 
# Thanks for reviewing this SciServer example notebook. You can use this notebook as a template to develop your own notebooks, but please do so in a copy rather than in the original example notebook.
# As you begin to use any of our SciServer modules in your own notebooks, consult the SciServer scripting documentation at http://www.sciserver.org/docs/sciscript-python/SciServer.html (link opens in a new window).
# 
# If you have questions, please email the SciServer helpdesk at sciserver-helpdesk@jhu.edu.

# In[ ]:




