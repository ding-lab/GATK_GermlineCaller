# Start local instance of cromwell server (rather than relying on default MGI production server)
# This gets around problems with database queries circa summer 2019
# With this server running, database queries are to http://localhost:8000
# Only one server should be running at once.  It should be run after `gsub` (after 0_start_docker.sh)
# and should exit when the docker container exits

CONFIG="/home/m.wyczalkowski/lib/cromwell-jar/server.cromwell.config"
CROMWELL="/usr/local/cromwell/cromwell-47.jar"

/usr/bin/java -Dconfig.file=$CONFIG -jar $CROMWELL server > /dev/null &

#/usr/bin/java -Dconfig.file=/home/m.wyczalkowski/lib/cromwell-jar/cromwell.truststore -jar /usr/local/cromwell/cromwell-47.jar server >/dev/null & 

