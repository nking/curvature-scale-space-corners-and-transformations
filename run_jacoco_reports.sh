#!/bin/bash

########################################################################
#  script to create coverage reports using jacoco/emma
#
########################################################################

ant clean compile compileTests

rm -rf bin/instr-classes
mkdir bin/instr-classes

cwd=`pwd`

jacocoJar=${JACOCO_HOME}/lib/jacocoant.jar

emmaJar=${cwd}/lib/emma/emma.jar

jarFiles3=$ANT_HOME/lib/ant.jar:~/.m2/repository/cglib/cglib-nodep/2.1_3/cglib-nodep-2.1_3.jar:$JAVA_HOME/lib/tools.jar:${cwd}/lib/junit.jar

classPath=$emmaJar:$jarFiles3

logLevel="INFO"
echo "making a logging.properties with level=${logLevel}"
rm bin/logging.properties

clsPath=.:${cwd}/bin/classes:bin/classes:$classPath

cd bin/test-classes

echo pwd >> ${cwd}/test_log.txt
echo >> ${cwd}/test_log.txt

filelist=`find . -name "*.class" -type f  | grep -v submodules | grep -v -F $ | grep Test| grep -v Aspect | grep -v TUtil | grep -v ForTests | grep -v dtrace | sed 's/\.class//g' | sed 's/\.\///g' | sed 's/\//\./g'`

rm ${cwd}/test_log.txt

echo "CLASSPATH=$clsPath" >> ${cwd}/test_log.txt
echo "" >> ${cwd}/test_log.txt
echo "FILELIST=$filelist" >> ${cwd}/test_log.txt
echo "" >> ${cwd}/test_log.txt

for file in $filelist; do
   echo "run test on $file" >> ${cwd}/test_log.txt

    java -javaagent:${JACOCO_HOME}/lib/jacocoagent.jar=destfile=${file}_jacoco.exec \
        -XX:-UseSplitVerifier -Xms512m -Xmx1024m -DisTesting=true \
        -Djava.util.logging.config.file=${cwd}/bin/logging.properties \
        -cp $clsPath org.junit.runner.JUnitCore ${file} >> ${cwd}/test_log.txt

done

cd ${cwd}

mkdir bin/jacoco

ant merge-exec-files -f ant_run_jacoco.xml

############################################################################
echo ""
echo "create the jacoco report"
############################################################################

ant report -f ant_run_jacoco.xml

