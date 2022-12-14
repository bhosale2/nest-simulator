#ifndef POLYMERINTEGRATOR_H_
#define POLYMERINTEGRATOR_H_

#include "Vector3.h"
#include "Matrix3.h"
#include "Rod.h"
#include "RodExternalForces.h"
#include "RodBoundaryConditions.h"
#include "Interaction.h"
//#include "MRAGProfiler.h"
#include "SpeedFunctions.h"
#include "ExternalContact.h"

// Symplectic integrator for a Polymer
class PolymerIntegrator
{
protected:
	typedef std::vector<Vector3> VV3;
	typedef std::vector<Matrix3> VM3;
	typedef std::vector<REAL> VREAL;
	typedef std::vector<bool> VBOOL;
	typedef std::vector<int> VINT;
	typedef std::vector<Rod*> Vrodptr;
	typedef std::vector<ExternalForces*> Vefptr;
	typedef std::vector<RodBC*> Vbcptr;
	typedef std::vector<Interaction*> Vinterptr;
	typedef std::vector<ExternalContact*> Vecptr;

	REAL time_old;

	Vrodptr rodptrs;
	Vefptr efptrs;
	Vbcptr bcptrs;
	Vinterptr interptrs;
	Vecptr ecptrs;

	void _computeForces(const REAL time, const int step);
	void _computeAccelerations(const REAL time, const int step);
	void _v_update_Q(const REAL coeffDt, Rod* rod);
	void _v_update_x(const REAL coeffDt, Rod* rod);
	void _v_update_v(const REAL coeffDt, Rod* rod);
	void _v_update_w(const REAL coeffDt, Rod* rod);

public:
	PolymerIntegrator(Vrodptr& _rodptrs, Vefptr& _efptrs, Vbcptr& _bcptrs, Vinterptr& _interptrs, Vecptr& _ecptrs) :
		time_old(0.0), rodptrs(_rodptrs), efptrs(_efptrs), bcptrs(_bcptrs), interptrs(_interptrs), ecptrs(_ecptrs) {}
	virtual ~PolymerIntegrator(){}

	inline Vrodptr& getRods(){ return rodptrs; }
	inline Vefptr& getExternalForces(){ return efptrs; }
	inline Vbcptr& getBoundaryConditions(){ return bcptrs; }
	inline Vinterptr& getInteractions(){ return interptrs; }
	inline Vecptr& getExternalContact(){ return ecptrs; }

	virtual REAL integrate(const REAL time, const REAL dt, const int step) = 0;
};

#endif


