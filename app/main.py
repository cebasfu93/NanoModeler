# -*- coding: utf-8 -*-

import pika
import random
import time
import os
import NanoModeler as nm
import json
import io
import base64
import sys
import traceback
import timeit
import logging

class LogStreamString(object):
	def __init__(self):
		self.logs = ''

	def write(self, data):
		self.logs += str(data)

	def flush(self):
		pass

	def clear(self):
		self.logs = ''

	def __str__(self):
		return self.logs.strip()

log_job_stream = LogStreamString()

def createFileHandler(content,mandatory=True):
	# check for empty content
	if not (isinstance(content, str) and (content and content.strip())):
		if mandatory:
			raise ValueError('Empty content not allowed')
		return None
	readFile = None
	# decode string from base64
	try:
		rawStringContent = base64.b64decode(content).decode('utf-8')
		readFile = io.StringIO(rawStringContent)
		readFile.name = 'temporary file'
	except Exception as e:
		print('createFileHandler failed: %s' % (str(e)))
		if mandatory:
			raise e
	finally:
		return readFile

def parseJson(content):
	resultJson = {} 
	if isinstance(content, str):
		return json.loads(content)
	elif isinstance(content, (bytes, bytearray)):
		return json.loads(content.decode('utf-8'))
	else:
		raise TypeError('content must be a string') 


def floatParameter(value,defaultValue):
	if(value is None):
		return float(defaultValue)
	try:
		return float(value)
	except TypeError:
		return defaultValue

def integerParameter(value,defaultValue):
	if(value is None):
		return int(defaultValue)
	try:
		return int(value)
	except TypeError:
		return defaultValue
	
def positiveIntegerParameter(value,defaultValue):
	if(value is None):
		return int(defaultValue)
	try:
		if int(value) < 0:
			return int(defaultValue)
		return int(value)
	except TypeError:
		return defaultValue	

def nvlParameter(value,defaultValue):
	if(value is None):
		return defaultValue
	return value

def parseParameters(request):
	parameters = {}

	inputParams	= request.get('parameters')

	# files
	parameters['LIG1_FILE'] = createFileHandler(inputParams.get('lig1'))
	parameters['LIG2_FILE'] = createFileHandler(inputParams.get('lig2'),False)
	parameters['CORE'] = createFileHandler(inputParams.get('core'))
	parameters['FRCMOD'] = createFileHandler(inputParams.get('frcMod'),False)

	# lists not mandatory
	parameters['CAP1'] = inputParams.get('cap1',[])
	parameters['CAP2'] = inputParams.get('cap2',[])

	# strings
	parameters['MORPHOLOGY'] = nvlParameter(inputParams.get('morphology'),'random')

	#integer positive
	parameters['STRIPES'] = positiveIntegerParameter(inputParams.get('stripes',1),1)
	parameters['RSEED'] = positiveIntegerParameter(inputParams.get('rseed',666),666)

	# float
	parameters['LIG1_FRAC'] = floatParameter(inputParams.get('lig1Frac',1.0),1.0)

	# integer
	parameters['LIG1_C'] = integerParameter(inputParams.get('lig1C',0),0)
	parameters['LIG1_S'] = integerParameter(inputParams.get('lig1S',0),0)
	parameters['LIG2_C'] = integerParameter(inputParams.get('lig2C',0),0)
	parameters['LIG2_S'] = integerParameter(inputParams.get('lig2S',0),0)

	#boolean
	parameters['ELONGATED'] = inputParams.get('elongated',False)

	return parameters

def processNanoModeler(request):
	parameters = parseParameters(request)
	return nm.NanoModeler(
		LIG1_FILE=parameters['LIG1_FILE'],
		CAP1=parameters['CAP1'],
		LIG1_C=parameters['LIG1_C'],
		LIG1_S=parameters['LIG1_S'],

		LIG1_FRAC=parameters['LIG1_FRAC'],
		MORPHOLOGY=parameters['MORPHOLOGY'],
		RSEED=parameters['RSEED'],
		STRIPES=parameters['STRIPES'],

		LIG2_FILE=parameters['LIG2_FILE'],
		CAP2=parameters['CAP2'],
		LIG2_C=parameters['LIG2_C'],
		LIG2_S=parameters['LIG2_S'],

		FRCMOD=parameters['FRCMOD'],

		CORE=parameters['CORE'],
		ELONGATED=parameters['ELONGATED'])

def on_request(ch, method, props, body):
	startTime = timeit.default_timer()
	uuid = props.correlation_id

	logging.info('start serving request - (%s)' % (uuid))

	time.sleep(5)

	resultJob = (0,None)

	errorMessage = ''

	# run job
	try:
		request = parseJson(body)
		resultJob = processNanoModeler(request)
	except Exception as e:
		print(traceback.format_exc())
		errorMessage = 'An error occurred while running job: %s' % (str(e))
		resultJob = (0,None)

	# read logs
	logs = str(log_job_stream)
	# reset logs
	log_job_stream.clear()

	version = nm.getVersion()
	status  = 'Completed' if resultJob[0] == 1 else 'Error'
	message = logs if resultJob[0] == 1 else '%s \nJob logs:\n%s' % (errorMessage,logs)

	resultObj = {}
	resultMessage = ''
	
	try:
		result  = base64.b64encode(resultJob[1]).decode('utf-8') if resultJob[0] == 1 else None

		resultObj['uuid'] = uuid
		resultObj['version'] = version
		resultObj['status'] = status
		resultObj['message'] = message
		resultObj['result'] = result
	except Exception as e:
		resultObj = {}
		resultObj['uuid'] = uuid
		resultObj['version'] = version
		resultObj['status'] = 'Error'
		resultObj['message'] = 'An error occurred while dumping result file: %s' % (str(e))
		resultObj['result'] = None
	finally:
		stopTime = timeit.default_timer()
		resultObj['execTimeMs'] = int(1000.0 * (stopTime - startTime))
		resultMessage = json.dumps(resultObj)

	logging.info('completed request - (%s) [%d ms]' % (uuid,resultObj['execTimeMs']))

	ch.basic_publish(exchange='',
				 routing_key=queue_response,
				 properties=pika.BasicProperties(correlation_id = uuid),
				 body=resultMessage);
	ch.basic_ack(delivery_tag = method.delivery_tag)

if __name__ == '__main__':

	nanoModelerLogger = logging.getLogger('nanomodeler')
	nanoModelerLogger.setLevel(logging.DEBUG)
	
	ch = logging.StreamHandler(log_job_stream)
	ch.setLevel(logging.INFO)

	### Optionally add a formatter
	formatter = logging.Formatter('%(levelname)s - %(message)s')
	ch.setFormatter(formatter)

	### Add the console handler to the logger
	nanoModelerLogger.addHandler(ch)
	
	logging.basicConfig(format= '%(message)s')

	logging.info('Starting pyconsumer daemon')

	username = os.getenv('RABBITMQ_USER', 'consumer')
	password = os.getenv('RABBITMQ_PASS', 'consumer')
	hostname = os.getenv('RABBITMQ_HOST', 'rabbitmq')

	queue_response = 'mdv.job.response'
	queue_listen = 'mdv.job.request'

	logging.info('Connection to rabbitmq service')

	credentials = pika.PlainCredentials(username, password)
	parameters = pika.ConnectionParameters(hostname, 5672, '/', credentials)

	connection = pika.BlockingConnection(parameters)

	logging.info('Connection enstablished with rabbitmq service')

	channel = connection.channel()

	channel.queue_declare(queue=queue_response,durable=True)
	channel.queue_declare(queue=queue_listen,durable=True)

	channel.basic_qos(prefetch_count=1)
	channel.basic_consume(on_request, queue=queue_listen)

	logging.info('Listening... Awaiting requests')

	channel.start_consuming()