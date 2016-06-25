/**
 \mainpage

 The ESP system support the use of machine learning for the analysis of
 real-time sensor data. ESP examples consist of code that specifies the
 application of machine learning to a particular application domain. An ESP
 example consists of multiple parts:

 \li Machine learning pipeline for a given application. This is an instance of
 GestureRecogntionPipeline from the GRT (http://nickgillian.com/grt/) specified
 using usePipeline().
 \li An input stream from which to receive live sensor data. This is an
 IStream instance specified by a call to useInputStream().
 \li Zero or one output streams to which to send the predictions generated by
 the machine learning pipeline. This is an OStream instance specified by
 a call to useOutputStream().
 \li Zero or more calibrators, which allow for the calibration of incoming
 sensor data based on the analysis of calibration samples collected by the
 user. Specified by the Calibrator instance passed to useCalibrator().
 \li Zero or more tuneable parameters, which allow the user to tune the value
 of various aspects of the application. Specified by the registerTuneable()
 functions in the tuneable.h file.
 \li Zero or one training samples checkers, to provide the user with feedback
 on the quality of their training data. Specified by a call to 
 useTrainingSampleChecker().
 \li Optional advice to the user on the collection of training data, using
 useTrainingDataAdvice().

 When compiled against the ESP framework and run, these examples generate a
 graphical interface that allows the user to visualize live sensor data,
 collect calibration samples, training, and test data, train the machine
 learning pipeline, adjust the tuneable parameters, and use the example in a
 larger interactive project.

 For an example, see: https://github.com/damellis/sensors/blob/master/Xcode/ESP/src/user_accelerometer_gestures.h
 */

#pragma once

#include "GRT/GRT.h"
#include "calibrator.h"
#include "iostream.h"
#include "tuneable.h"
#include "training.h"

using namespace GRT;
