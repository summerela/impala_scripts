##install R packages
install.packages('rJava')
install.packages('RImpala')

library(RImpala)
library(RCurl)

##make sure java is ready for this
system(paste("sudo R CMD javareconf"), wait=TRUE)

##install impala jdbc jars by downloading and unzipping
download.file(url="https://downloads.cloudera.com/impala-jdbc/impala-jdbc-0.5-2.zip", destfile="/Users//summerrae/Documents/", 
              method="curl")
system(paste("export CLASSPATH=/Users/summerrae/Documents/impala-jdbc-0.5-2/:$CLASSPATH"))

#initialize impala jdbc jars
rimpala.init(libs="/Users/summerrae/Documents/impala-jdbc-0.5-2")

#establish connection
rimpala.connect("glados19", "21050")

##see what's in here
rimpala.showdatabases()

##select a database to use
rimpala.usedatabase("public_hg19") 

##what tables are in here? 
rimpala.showtables()

##lets peek at clinvar
rimpala.describe("clinvar") 

test_query = rimpala.query("Select * from kaviar limit 3")

rimpala.close()

