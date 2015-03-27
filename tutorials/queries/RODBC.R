install.packages("RODBC")
library(RODBC)

#connect using the DSN you created above
conn <- odbcConnect("Impala DSN")

#view available tables
sqlTables(conn)

##read a table into a data frame
crimedat <- sqlFetch(myconn, "Crime")


#close the connection
close(conn)
