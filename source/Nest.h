/*
 * Nest.h
 *
 *  Created on: Jun 22, 2014
 *      Author: mgazzola
 */

#ifndef NEST_H_
#define NEST_H_

#include "UsualHeaders.h"
#include "Test.h"
#include "ArgumentParser.h"
#include "RodInitialConfigurations.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"

using namespace std;

class Nest: public Test
{
protected:
        int NOR;        

public:

	Nest(const int argc, const char ** argv);
	~Nest(){};

	void run();
	void paint(){};
};



#endif /* NEST_H_ */
