/*
 * Collision_test.h
 *
 *  Created on: Jun 22, 2014
 *      Author: mgazzola
 */

#ifndef COLL_H_
#define COLL_H_

#include "UsualHeaders.h"
#include "Test.h"
#include "ArgumentParser.h"
#include "RodInitialConfigurations.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"

using namespace std;

class Coll_test: public Test
{
protected:
        int NOR;        

public:

	Coll_test(const int argc, const char ** argv);
	~Coll_test(){};

	void run();
	void paint(){};
};



#endif /* Coll_test_H_ */
