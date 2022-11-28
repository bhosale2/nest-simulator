/*
 * Mindlin_friction.h
 *
 *  Created on: Jun 22, 2014
 *      Author: mgazzola
 */

#ifndef MIND_H_
#define MIND_H_

#include "UsualHeaders.h"
#include "Test.h"
#include "ArgumentParser.h"
#include "RodInitialConfigurations.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"

using namespace std;

class Mindlin_friction: public Test
{
protected:
        int NOR;        

public:

	Mindlin_friction(const int argc, const char ** argv);
	~Mindlin_friction(){};

	void run();
	void paint(){};
};



#endif /* NEST_H_ */
