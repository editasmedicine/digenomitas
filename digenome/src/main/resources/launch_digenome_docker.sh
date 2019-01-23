#!/bin/bash

# bash script to run digenome docker containers on local or batch
# provide:
#   -r req
#   -b analysis_sheet_bucket
#   -k analysis_sheet_key
# examples:
# ./launch_digenome_docker.sh local -r ngs-req-640 -b editas-fulcrumgeonomics -k digenome-test-data/run-folder/sample-sheet.xlsx -d -f ../../../../pipeline/test_data/event.json

EVENT_TEMPLATE_FILE="$(dirname $0)/digenome-1.0.json"
TEST_TEMPLATE_FILE="$(dirname $0)/../../../../pipeline/test_data/event.json"
DOCKER_IMAGE="digenome"

# for SQS executions
SQS_QUEUE_NAME="analysis-pipeline"
SQS_URL_PREFIX="https://queue.amazonaws.com/422765174246/"

BATCH_TOPIC="arn:aws:sns:us-east-1:422765174246:sequencingAnalysisStatus"
SNS_URL_PREFIX="https://console.aws.amazon.com/sns/v2/home?region=us-east-1#/topics/"

SYNOPSIS_STRING="usage: ./launch_digenome_docker.sh [ local | batch | sqs] -r <req> -b <s3_bucket> -k <analysis_sheet_key> [options]"
HELPSTRING=\
"Digenome Run Script
  SYNOPSIS
    $SYNOPSIS_STRING
  
  ARGUMENTS
    -r : req (string)
    -b : analysis sheet S3 bucket (string)
    -k : analysis sheet S3 key (string)
  
  OPTIONS
    -j : jira (string)
      jira issue for updates
    
    -n : ncpu (int)
      override number of vcpus for batch execution
      
    -m : memory (int)
      override memory for batch execution

    -f: file for json template (string)
      path to json file to override tempalte. Takes $EVENT_TEMPLATE_FILE
      
    -d: dev (boolean flag)
      runs dev version of the batch pipeline and send to .dev s3 locations
      
    -i: inspect (boolean flag)
      For local runs. Does not run the start CMD. Instead runs the container with an interactive bash session."

if [[ $1 == 'help' ]]; then
  printf "$HELPSTRING"
  exit 0
fi

# if "test" argument passed, run with default testing parameters
if [[ $2 == "test" ]]; then
  echo "running as test with defaul parameters"
  ./$0 $1 -r ngs-req-640-subsample -b editas-fulcrumgeonomics -k digenome-test-data/run-folder/sample-sheet.xlsx -f $TEST_TEMPLATE_FILE
  exit 0
fi

# Determined whether this is local or batch
runType=$1
if [[ $runType == "local" || $runType == "batch" ]]; then
  # continue
  shift
else
  echo $SYNOPSIS_STRING
  echo "use batch or local"
  exit 1
fi


req=""
bucket=""
key=""
jira=""
ncpu=""
memory=""
samples="["
dev=0

while getopts "j:n:m:r::b:k:f:d-:s:" opt; do
  case $opt in
    r)
      req=$OPTARG >&2;;
    b)
      bucket=$OPTARG >&2;;
    k)
      key=$OPTARG >&2;;
    s)
      samples+="$prefix$OPTARG"
      prefix=", ";;
    j)
      echo "adding jira issue: $OPTARG to event" >&2
      jira=$OPTARG;;
    n)
      echo "modifying ncpus to $OPTARG" >&2
      ncpu=$OPTARG;;
    m) 
      echo "modifying memory requirement to $OPTARG" >&2
      memory=$OPTARG;;
    f)
      echo "using your template file instead: $OPTARG" >&2
      EVENT_TEMPLATE_FILE=$OPTARG;;
    d) dev=1;;
    i) interactiveSession="/bin/bash";;
    \?)
      echo "Invalid option" >&2
      exit 1;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1;;
  esac
done
shift $((OPTIND -1))

samples=$samples"]"

if [[ $req == "" || $bucket == "" || $key == "" ]]; then
    echo "Must provide REQ, BUCKET, and KEY with appropriate flags: -r, -b, -k"
    exit 1
fi

# add req, bucket, and key
event=$(jq ".req=\"$req\"" $EVENT_TEMPLATE_FILE | jq ".analysisSheetLocation.bucket=\"$bucket\"" | jq ".analysisSheetLocation.key=\"$key\"")

# add samples if provided
if [[ $samples != "[]" ]]; then
  event=$(echo $event | jq ".analysisParameters.process_samples=$samples")
fi

# add options
if [[ $jira != "" ]]; then
  event=$(echo $event | jq ".jira.issue=\"$jira\"")
fi
if [[ $ncpu != "" ]]; then
  event=$(echo $event | jq ".batchJob.containerOverrides.vcpus=$ncpu")
fi
if [[ $memory != "" ]]; then
  event=$(echo $event | jq ".batchJob.containerOverrides.memory=$memory")
fi

if [[ $dev == 1 ]]; then
    echo "executing as dev version"
    resultsBucket=$(echo $event | jq -r ".resultsData.bucket")
    parquetBucket=$(echo $event | jq -r ".parquetData.bucket")
    rawDataBucket=$(echo $event | jq -r ".rawData.bucket")
    devExtension="-dev"
    devBucketExtension=".dev"
    event=$(echo $event | jq ".resultsData.bucket=\"$resultsBucket$devBucketExtension\"" | jq ".parquetData.bucket=\"$parquetBucket$devBucketExtension\"")
    if [[ ${rawDataBucket} != *"$devBucketExtension"* ]]; then
      event=$(echo $event | jq ".rawData.bucket=\"$rawDataBucket$devBucketExtension\"")
    fi
fi

echo "Executing on $runType environment"
echo "$event" | jq '.' #jq used for coloring
echo "continue? [y/N]:"
read response

if [[ $response == "y" || $response == "Y" ]]; then
  echo "executing"
  if [[ $runType = "local" ]]; then
    docker run --rm -ti -v /home/ec2-user/efs/reference:/root/reference \
                      -e EVENT="$event" \
                      $DOCKER_IMAGE $interactiveSession
  elif [[ $runType = "batch" ]]; then
    response=$(aws stepfunctions start-execution \
      --state-machine-arn  "arn:aws:states:us-east-1:422765174246:stateMachine:batchSequencing"$devExtension \
      --input "$event")
    echo $response | jq '.'
    executionArn=$(echo $response | jq -r '.executionArn')
    SFN_URL="https://console.aws.amazon.com/states/home?region=us-east-1#/executions/details/"
    echo "See SFN execution at:"
    echo $SFN_URL$executionArn
  else
    # post to sqs
    echo 'posting to sqs'
    sqs_url=$SQS_URL_PREFIX$SQS_QUEUE_NAME$devExtension
    echo $sqs_url
    echo $(aws sqs send-message --queue-url $sqs_url --message-body "$event") | jq '.'
    echo "subscribe to SNS topic $BATCH_TOPIC$devExtension at $SNS_URL_PREFIX$BATCH_TOPIC$devExtension for updates"
    exit 0
  fi
else
  echo "aborting."
  exit 0
fi
