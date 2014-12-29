//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2010-2012 Alessandro Tasora
// Copyright (c) 2013 Project Chrono
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be 
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

///////////////////////////////////////////////////
//      
//   Demo code about   
//   
//     - collisions and contacts 
//
//       (This is just a possible method of integration
//       of Chrono::Engine + Irrlicht: many others 
//       are possible.)
//     
//	 CHRONO   
//   ------
//   Multibody dinamics engine 
//   
// ------------------------------------------------ 
//             www.deltaknowledge.com
// ------------------------------------------------ 

//#include "physics/ChBodyEasy.h"
#include "physics/ChContactContainer.h"
#include "collision/ChCModelBulletBody.h"
//#include "core/ChTimer.h"
//#include "core/ChRealtimeStep.h"
//#include "assets/ChTexture.h"
#include "unit_IRRLICHT/ChIrrApp.h"
#include <cstring>
#include <fstream>
#include <sstream>
#include <time.h>
#include <cstdlib>
//#include <map>

//*************** chrono parallel
#include <stdio.h>
#include <vector>
#include <cmath>

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

#include "chrono_utils/ChUtilsCreators.h"  //Arman: why is this
#include "chrono_utils/ChUtilsInputOutput.h" //Arman: Why is this
#include "chrono_utils/ChUtilsGenerators.h"

#ifdef CHRONO_PARALLEL_HAS_OPENGL2
#include "chrono_opengl/ChOpenGLWindow.h"
#endif
//***********************************
// Use the namespace of Chrono

using namespace chrono;
using namespace chrono::collision;

// Use the main namespaces of Irrlicht
using namespace irr;
using namespace core;
using namespace scene;
using namespace video;
using namespace io;
using namespace gui;
using namespace std;

enum DriveType {
	ACTUATOR,
	KINEMATIC
};
DriveType driveType;
//******************* Initialize attributes *******************
const double rhoF = 1000;
const double rhoR = 917;
const double rhoPlate = 1000;
const double porosity = .35; //set value of porosity
const double minPorosity = 0.2595;
const double mu_Viscosity = .001;//.1;
const ChVector<> surfaceLoc = ChVector<>(0, .04, -.08);

double mradius = 0.4;
double collisionEnvelop = .04 * mradius;
int numLayers = 3;
const double shipVelocity = 5.4;//.27;//1; //arman modify
const double timePause = 1;//1;//0.2; //arman modify : Time pause != 0 causes the actuator to explode
const double timeMove = 2.5;

// ** box and ship locations **
const double ship_width = 4;
const double box_X = ship_width, box_Y = 10, box_Z = .4;
ChVector<> shipInitialPos;
double shipInitialPosZ = 0;
ChVector<> hdim = ChVector<>(16, 12, 21.2); //domain dimension
ChVector<> boxMin = ChVector<>(-.8, -6, -2.4); //component y is not really important
//**************************************************************
ChSharedBodyPtr shipPtr;
ChSharedBodyPtr bin;
//**************************************************************
void MySeed(double s = time(NULL)) {
	 srand(s);
}
double MyRand() {
	return float(rand()) / RAND_MAX;
}
void Calc_Hydrodynamics_Forces(ChVector<> & F_Hydro, ChVector<> & forceLoc, ChVector<> & T_Drag,
		ChBody* mrigidBody, ChSystemParallel& mphysicalSystem, const chrono::ChVector<>& freeSurfaceLocation) {
	F_Hydro = ChVector<>(0,0,0);
	forceLoc = ChVector<>(0,0,0);
	T_Drag = ChVector<>(0,0,0);


	// ***** calculation of force
	ChVector<> freeSurfaceNormal = -mphysicalSystem.Get_G_acc(); // The default acceleration is (0, -9.81, 0);
	double g = freeSurfaceNormal.Length();
	freeSurfaceNormal.Normalize();

	ChVector<> bodyCtr = mrigidBody->GetPos();
	ChVector<> dist3 = bodyCtr - freeSurfaceLocation;
	double dist = dist3.Dot(freeSurfaceNormal); // distance of the sphere center from the fluid surface

	//mrigidBody->GetCollisionModel()->GetSafeMargin(); //this only works for sphere
	double rad = mradius;
	//****************** Buoyancy Force
	ChVector<> F_Buoyancy = ChVector<>(0,0,0);
	forceLoc = bodyCtr;
	if (dist < -rad) {
		double V = 4.0 / 3 * CH_C_PI * pow(rad, 3);
		F_Buoyancy = V * rhoF * g * freeSurfaceNormal;
		forceLoc = bodyCtr;
	} else if (dist < rad) {
		double h = rad - dist;
		double V = CH_C_PI * h * h / 3 * (3 * rad - h);
		F_Buoyancy = V * rhoF * g * freeSurfaceNormal;
		double distFromCenter = 3.0 / 4 * pow(2 * rad - h, 2) / (3 * rad - h); 	// http://mathworld.wolfram.com/SphericalCap.html -->
																				// Harris and Stocker 1998, p. 107 (Harris, J. W. and Stocker,
																				// H. "Spherical Segment (Spherical Cap)." §4.8.4 in Handbook of
																				// Mathematics and Computational Science. New York: Springer-Verlag, p. 107, 1998.)
		forceLoc = bodyCtr + distFromCenter * (-freeSurfaceNormal);
	}
	// "dist > rad" --> outside of water
	//****************** Drag Force and Torque
	double Cd = 0.4;
	ChVector<> vel = mrigidBody->GetPos_dt();
	ChVector<> F_Drag = ChVector<>(0,0,0);
//	if (dist < rad) {
//		double A_ref = 0.5 * CH_C_PI * rad * (rad - dist);
//		double multDrag = 1;
//		if (mphysicalSystem.GetChTime() < timePause) {
//			multDrag = 1;
//		} else {
//			multDrag = 1;
//		}
//		F_Drag = multDrag * (-6.0 * CH_C_PI * mu_Viscosity * rad * vel
//					-0.5 * rhoF * Cd * vel.Length() * vel);
//		T_Drag = -8.0 * CH_C_PI * mu_Viscosity * pow(rad, 3) * mrigidBody->GetWvel_par(); // in parent, i.e. absoute, reference frame.
//	}
	//****************** Total Force
	F_Hydro = F_Buoyancy + F_Drag; // it is assumed that F_Drag is applied at the buoyancy center
}
//**********************************
void create_hydronynamic_force(ChBody* mrigidBody, ChSystemParallel& mphysicalSystem, const chrono::ChVector<>& freeSurfaceLocation, bool createForce) {
	// ***** insertion of force
	ChSharedPtr<ChForce> hydroForce;
	ChSharedPtr<ChForce> hydroTorque;

//	string forceTag("hydrodynamics_force");
	char forceTag[] = "hydrodynamics_force";
	char torqueTag[] = "hydrodynamics_torque";
	hydroForce = mrigidBody->SearchForce(forceTag);
	hydroTorque = mrigidBody->SearchForce(torqueTag);

	//********** create force if needed **********
	if (createForce) {
		if (hydroForce.IsNull()) {
			hydroForce = ChSharedPtr<ChForce>(new ChForce);
			hydroForce->SetMode(FTYPE_FORCE); // no need for this. It is the default option.
			mrigidBody->AddForce(hydroForce);
			// ** or: hydroForce = ChSharedPtr<ChForce>(new ChForce());
			hydroForce->SetName(forceTag);
		}
		if (hydroTorque.IsNull()) {
			hydroTorque = ChSharedPtr<ChForce>(new ChForce);
			hydroTorque->SetMode(FTYPE_TORQUE);
			mrigidBody->AddForce(hydroTorque);
			// ** or: hydroForce = ChSharedPtr<ChForce>(new ChForce());
			hydroTorque->SetName(torqueTag);
		}
	}
	//********** update force magnitude **********
	if (!hydroForce.IsNull() || !hydroTorque.IsNull()) {
		ChVector<> F_Hydro;
		ChVector<> forceLoc;
		ChVector<> T_Drag;

		Calc_Hydrodynamics_Forces(F_Hydro, forceLoc, T_Drag, mrigidBody, mphysicalSystem, freeSurfaceLocation);

		hydroForce->SetVpoint(forceLoc);
		hydroForce->SetMforce(F_Hydro.Length());
		F_Hydro.Normalize();
		hydroForce->SetDir(F_Hydro);

		hydroTorque->SetMforce(T_Drag.Length());
		T_Drag.Normalize();
		hydroTorque->SetDir(T_Drag);
	}
}
//**********************************
void calc_ship_contact_forces(ChSystemParallelDVI& mphysicalSystem, ChVector<> & mForce, ChVector<> & mTorque) {
	mForce = ChVector<>(0,0,0);
	mTorque = ChVector<>(0,0,0);
	real3 myForce;
	ChContactContainer* container  = (ChContactContainer *) mphysicalSystem.GetContactContainer();
//	map<ChBody*, ChVector<> > m_forces;
//	map<ChBody*, ChVector<> > m_torques;
//	ChVector<> mForce;
//	ChVector<> mTorque;

	unsigned int bodyID = shipPtr->GetId();// (((ChBody)shipPtr)->GetId());
	int offset = 3; //Arman modify this
	if (mphysicalSystem.GetSettings()->solver.solver_mode != SLIDING) {
		offset = 3;
	} else if (mphysicalSystem.GetSettings()->solver.solver_mode != SPINNING) {
		offset = 6;
	}

	real3 gam(0,0,0);
	real3 gammaRoll(0,0,0);
	double dT = mphysicalSystem.GetStep();

	for (int i = 0; i < mphysicalSystem.data_manager->num_contacts; i++) {
		int2 ids = mphysicalSystem.data_manager->host_data.bids_rigid_rigid[i];
		real3 U = mphysicalSystem.data_manager->host_data.norm_rigid_rigid[i];
		if (ids.x != bodyID && ids.y != bodyID) continue;
		if (mphysicalSystem.GetSettings()->solver.solver_mode == NORMAL) {
			gam.x = mphysicalSystem.data_manager->host_data.gamma_data[i * offset + 0];
			myForce += (-gam.x * dT * U);
		} else if (mphysicalSystem.GetSettings()->solver.solver_mode == SLIDING) { //Arman check this
			gam.x = mphysicalSystem.data_manager->host_data.gamma_data[i * offset + 0];
			gam.y = mphysicalSystem.data_manager->host_data.gamma_data[i * offset + 1];
			gam.z = mphysicalSystem.data_manager->host_data.gamma_data[i * offset + 2];
			real3 V, W;
			Orthogonalize(U, V, W);
			myForce += dT * (-U * gam.x - V * gam.y - W * gam.z);
//			//*** torque related
//		    bool2 active = contact_active[index];
//		    if (active.x != 0) {
//		       real3 T3, T4, T5;
//		       Compute_Jacobian(rot[index], U, V, W, ptA[index], T3, T4, T5);
//		       updateO[index] = T3 * gam.x + T4 * gam.y + T5 * gam.z;
//		    }
//		    if (active.y != 0) {
//		       real3 T6, T7, T8;
//		       Compute_Jacobian(rot[index + num_contacts], U, V, W, ptB[index], T6, T7, T8);
//		       //updateV[index + num_contacts] = U * gam.x + V * gam.y + W * gam.z;
//		       updateO[index + num_contacts] = -T6 * gam.x - T7 * gam.y - T8 * gam.z;
//		    }
		} else if (mphysicalSystem.GetSettings()->solver.solver_mode == SPINNING) { //Arman check this
			gam.x = mphysicalSystem.data_manager->host_data.gamma_data[i * offset + 0];
			gam.y = mphysicalSystem.data_manager->host_data.gamma_data[i * offset + 1];
			gam.z = mphysicalSystem.data_manager->host_data.gamma_data[i * offset + 2];
			real3 V, W;
			Orthogonalize(U, V, W);
			myForce += (-U * gam.x - V * gam.y - W * gam.z) * dT;
//			//*** torque related
//			if (active.x != 0) {
//			 real4 quat = rot[index];
//			 real3 T3, T4, T5, TA, TB, TC;
//			 Compute_Jacobian(quat, U, V, W, ptA[index], T3, T4, T5);
//			 Compute_Jacobian_Rolling(quat, U, V, W, TA, TB, TC);
//
//			 updateO[index] = T3 * gam.x + T4 * gam.y + T5 * gam.z - TA * gam_roll.x - TB * gam_roll.y - TC * gam_roll.z;
//			}
//			if (active.y != 0) {
//			 real4 quat = rot[index + num_contacts];
//			 real3 T6, T7, T8, TA, TB, TC;
//			 Compute_Jacobian(quat, U, V, W, ptB[index], T6, T7, T8);
//			 Compute_Jacobian_Rolling(quat, U, V, W, TA, TB, TC);
//			 //updateV[index + num_contacts] = U * gam.x + V * gam.y + W * gam.z;
//			 updateO[index + num_contacts] = -T6 * gam.x - T7 * gam.y - T8 * gam.z + TA * gam_roll.x + TB * gam_roll.y + TC * gam_roll.z;
//			}
		}
	}


	//extra torque term
//	real3 bodyCenter = 	mphysicalSystem.data_manager->host_data.ObA_rigid[bodyID];
//	real3 contactPt = 	mphysicalSystem.data_manager->host_data.cpta_rigid_rigid[bodyID];
//	mTorque = ChVector<>(contactPt.x - bodyCenter.x, contactPt.y - bodyCenter.y, contactPt.z - bodyCenter.z) % mForce;
	mForce = ChVector<>(myForce.x, myForce.y, myForce.z);
}
//***********************************
void CreateSphere(ChSystemParallelDVI& mphysicalSystem, ChVector<> pos) {
	static int sphereID = 0;

	double mmass = (4./3.)*CH_C_PI*pow(mradius,3)*rhoR;
	double minert = (2./5.)* mmass * pow(mradius,2);
	ChSharedBodyPtr mrigidBody(new ChBody(new ChCollisionModelParallel));

	// set moment of inertia (more realistic than default 1,1,1).
	mrigidBody->SetIdentifier(sphereID++);
	mrigidBody->SetMass(mmass);
	mrigidBody->SetInertiaXX(minert * ChVector<>(1, 1, 1));
	mrigidBody->SetPos(pos);
	mrigidBody->SetPos_dt(ChVector<>(0,0,0));
	mrigidBody->SetRot(ChQuaternion<>(1, 0, 0, 0));
	mrigidBody->SetBodyFixed(false);
	mrigidBody->SetCollide(true);


	ChSharedPtr<ChMaterialSurface> ballMat(new ChMaterialSurface);
	ballMat->SetFriction(0.4f);
	ballMat->SetCompliance(0.0);
	ballMat->SetComplianceT(0.0);
	ballMat->SetDampingF(0.2);
	mrigidBody->SetMaterialSurface(ballMat);

	//		mrigidBody->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelop); //envelop is .03 by default
	mrigidBody->GetCollisionModel()->ClearModel();
	utils::AddSphereGeometry(mrigidBody.get_ptr(), mradius);
	mrigidBody->GetCollisionModel()->BuildModel();

//	// optional, attach a texture for better visualization  //Arman irrlicht
//	ChSharedPtr<ChTexture> mtextureball(new ChTexture());
//	mtextureball->SetTextureFilename(GetChronoDataFile("../data/bluwhite.png"));
//	mrigidBody->AddAsset(mtextureball);

	mphysicalSystem.Add(mrigidBody);
}
//***********************************
void GenerateIceLayers_Rectangular(
		ChSystemParallelDVI& mphysicalSystem,
		ChVector<> boxMin,
		ChVector<> boxMax,
		double global_x, //global offset in x
		double global_y, //global offset in y
		double global_z, //global offset in z
		double expandR)
{
	int numColX = (boxMax.x - boxMin.x - expandR) / (2 * expandR);
	int numColZ = (boxMax.z - boxMin.z - expandR) / (2 * expandR);
	double spacing = 2 * expandR;
	for (int j = 0; j < numLayers; j++) {
		for (int i = 0; i < numColX; i++) {
			for (int k = 0; k < numColZ; k++) {
				// Create a ball that will collide with wall
				ChVector<> pos = ChVector<>(
						ChVector<>(global_x, global_y, global_z)
						+ ChVector<>((i+0.5) * spacing, j * spacing, (k+0.5) * spacing)
						+ .5 * (spacing - 2 * mradius) * ChVector<>(ChRandom(), ChRandom(), ChRandom())
						);
				CreateSphere(mphysicalSystem, pos);
			}
		}

	}
}

//***********************************
void GenerateIceLayers_Hexagonal(
		ChVector<> boxMin,
		ChVector<> boxMax,
		double global_x, //global offset in x
		double global_y, //global offset in y
		double global_z,
		double expandR,

		vector<double> & px, vector<double> & py, vector<double> & pz) //global offset in z
{
	int numColX = (boxMax.x - boxMin.x - 2 * expandR) / (2 * expandR);
	int numColZ = (boxMax.z - boxMin.z - 2 * expandR) / (sqrt(3.0) * expandR);
    double offset_x = 0, offset_z = 0;
    for (int j = 0; j < numLayers; j++) {
      double y = j * (sqrt(3.0) * expandR);
      //need to offset each alternate layer by radius in both x and z direction
      offset_x = offset_z = (j % 2 != 0) ? expandR : 0;
      double offset = 0;
      double x = 0, z = 0;
      for (int i = 0; i < numColX; i++) {
        for (int k = 0; k < numColZ; k++) {
			offset = (k % 2 != 0) ? expandR : 0;
			x = i * 2 * expandR + offset + (expandR + global_x + offset_x);
			z = k * (sqrt(3.0) * expandR)  + (expandR + global_z + offset_z);
			ChSharedBodyPtr mrigidBody(new ChBody(new ChCollisionModelParallel));
			ChVector<> pos = ChVector<>(x, y, z);
			px.push_back(pos.x);
			py.push_back(pos.y);
			pz.push_back(pos.z);
        }
      }
    }
}

void AddParticlesToSys(
		ChSystemParallelDVI& mphysicalSystem, vector<double> & px, vector<double> & py, vector<double> & pz) //global offset in z
{
	int numPart = px.size();
	for (int i = 0; i < numPart; i++) {
		ChVector<> pos = ChVector<>(px[i], py[i], pz[i]);
		ChSharedBodyPtr mrigidBody(new ChBody(new ChCollisionModelParallel));
		CreateSphere(mphysicalSystem, pos);
	}
}

// =============================================================================
// Create the Brash Ice
//
// Brash Ice consisting of identical spheres with specified radius and
// material properties; the spheres are generated in a number of vertical
// layers with locations within each layer obtained using Poisson Disk sampling,
// or HCP Packing, or Grid based initializaiton, thus ensuring that no two spheres
// are closer than twice the radius.
// =============================================================================
int CreateIceParticles(ChSystemParallel& mphysicalSystem)
{
	// -------------------------------------------
	// Create a material for the granular material
	// -------------------------------------------
	ChSharedPtr<ChMaterialSurface> mat_g(new ChMaterialSurface);

	mat_g->SetFriction(0.4f);
	mat_g->SetCompliance(0.0);
	mat_g->SetComplianceT(0.0);
	mat_g->SetDampingF(0.2);

	// ---------------------------------------------
	// Create a mixture entirely made out of spheres
	// ---------------------------------------------

	// Create the particle generator with a mixture of 100% spheres
	utils::Generator gen(&mphysicalSystem);
	utils::MixtureIngredientPtr& m1 = gen.AddMixtureIngredient(utils::SPHERE, 1.0);
	m1->setDefaultMaterialDVI(mat_g);
	m1->setDefaultDensity(rhoR);
	m1->setDefaultSize(mradius);

	// Ensure that all generated particle bodies will have positive IDs.
	int Id_g = 1;
	gen.setBodyIdentifier(Id_g);

	// ----------------------
	// Generate the particles
	// ----------------------

	double expandR = 1.05 * mradius;
	double iceThickness = numLayers * expandR * 2;
	double buttomLayerDY = rhoR / rhoF *  iceThickness;
	ChVector<> boxMinGranular = ChVector<>(boxMin.x, surfaceLoc.y - buttomLayerDY, boxMin.z);
//	ChVector<> hdimGranularHalf = 0.5 * hdim - ChVector<>(expandR);
	ChVector<> hdimGranularHalf = 0.5 * ChVector<>(hdim.x, iceThickness, hdim.z) - ChVector<>(expandR);
	ChVector<> centerGranular = boxMinGranular + (hdimGranularHalf + ChVector<>(expandR));

	printf("************************** Generate Ice, ButtomLayer_Y %f\n", buttomLayerDY);

	gen.createObjectsBox(utils::HCP_PACK, 2 * expandR, centerGranular, hdimGranularHalf); //REGULAR_GRID : HCP_PACK

	// Return the number of generated particles.
	return gen.getTotalNumBodies();
}

// =============================================================================
// Add customized attributes the set of bodies in the physical system that are
// in the range [idx_i, idx_j]
// =============================================================================
void AddCustomAttribute(ChSystemParallel& mphysicalSystem, int idx_i, int idx_j) {
	std::vector<ChBody*>::iterator myIter = mphysicalSystem.Get_bodylist()->begin() + idx_i;
	for (int i = idx_i; i < idx_j; i++) {
		create_hydronynamic_force(*myIter, mphysicalSystem, surfaceLoc, true);  //Arman : hydrodynamic forces
		myIter++;
	}
}
//***********************************
void create_system_particles(ChSystemParallelDVI& mphysicalSystem)
{
	//**************** sphere prob
	double expandR = mradius*1.05;

	double iceThickness = numLayers * expandR * 2;
	double buttomLayerDY = rhoR / rhoF *  iceThickness - mradius;

	ChVector<> center = 0.5 * hdim + boxMin;
	ChVector<> boxMax = boxMin + hdim;
	printf("************************** Generate Ice, ButtomLayer_Y %f\n", buttomLayerDY);
	double global_x = boxMin.x;
	double global_y = surfaceLoc.y - buttomLayerDY;
	double global_z = boxMin.z;

	int idxI = mphysicalSystem.Get_bodylist()->size();

//	GenerateIceLayers_Rectangular(mphysicalSystem,
//			boxMin, boxMax,
//			global_x, global_y, global_z,
//			expandR);

//	vector<double> px, py, pz;
//	GenerateIceLayers_Hexagonal(boxMin, boxMax, global_x, global_y, global_z, expandR, px, py, pz);
//	AddParticlesToSys(mphysicalSystem, px, py, pz);

	(void)CreateIceParticles(mphysicalSystem);


	int idxJ = mphysicalSystem.Get_bodylist()->size();
	AddCustomAttribute(mphysicalSystem, idxI, idxJ);

	//**************** bin and ship
	// IDs for the two bodies
	int binId = -200;
	int shipId = -201;

	// Create a common material
	ChSharedPtr<ChMaterialSurface> mat(new ChMaterialSurface);
	mat->SetFriction(0.4f);
	mat->SetDampingF(0.2f);

	// Create the containing bin (2 x 2 x 1)
	double hthick = .1;
	double hole_width = 1.2 * ship_width;
	double small_wall_Length = 0.5 * (hdim.x - hole_width);

	bin = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	bin->SetMaterialSurface(mat);
	bin->SetIdentifier(binId);
	bin->SetMass(1);
	bin->SetPos(ChVector<>(center.x, center.y, center.z));
	bin->SetRot(ChQuaternion<>(1, 0, 0, 0));
	bin->SetCollide(true);
	bin->SetBodyFixed(true);
	bin->GetCollisionModel()->ClearModel();

	utils::AddBoxGeometry(bin.get_ptr(), 0.5 * ChVector<>(hdim.x, hdim.y, hthick), ChVector<>(0, 0, 0.5 * hdim.z + 0.5*hthick));	//end wall
	utils::AddBoxGeometry(bin.get_ptr(), 0.5 * ChVector<>(hthick, hdim.y, hdim.z), ChVector<>(-0.5 * hdim.x - 0.5 * hthick, 0, 0));		//side wall
	utils::AddBoxGeometry(bin.get_ptr(), 0.5 * ChVector<>(hthick, hdim.y, hdim.z), ChVector<>(0.5 * hdim.x + 0.5 * hthick, 0, 0));	//side wall
	utils::AddBoxGeometry(bin.get_ptr(), 0.5 * ChVector<>(small_wall_Length, hdim.y, hthick), ChVector<>(-0.5 * hdim.x + 0.5*small_wall_Length, 0, -0.5 * hdim.z - 0.5*hthick)); 	//beginning wall 1
	utils::AddBoxGeometry(bin.get_ptr(), 0.5 * ChVector<>(small_wall_Length, hdim.y, hthick), ChVector<>(0.5 * hdim.x - 0.5*small_wall_Length, 0, -0.5 * hdim.z - 0.5*hthick)); //beginning wall 2
	utils::AddBoxGeometry(bin.get_ptr(), 0.5 * ChVector<>(7 * hdim.x, hthick, 7 * hdim.x), ChVector<>(0,-10,0)); //bottom bed
	bin->GetCollisionModel()->BuildModel();

	mphysicalSystem.AddBody(bin);

	//**************** create ship
	double boxMass = rhoPlate * box_X * box_Y * box_Z;
	double bI1 = 1.0 / 12 * boxMass * (pow(box_X, 2) + pow(box_Y, 2));
	double bI2 = 1.0 / 12 * boxMass * (pow(box_Y, 2) + pow(box_Z, 2));
	double bI3 = 1.0 / 12 * boxMass * (pow(box_X, 2) + pow(box_Z, 2));
	printf("mass %f I1 I2 I3 %f %f %f\n", boxMass, bI1, bI2, bI3);

	shipInitialPosZ = boxMin.z - .5 * box_Z;

	shipPtr = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	shipInitialPos = ChVector<>(center.x,  center.y, shipInitialPosZ);
	shipPtr->SetPos(shipInitialPos);
	shipPtr->SetRot(ChQuaternion<>(1,0,0,0));
	shipPtr->SetMaterialSurface(mat);
	shipPtr->SetPos_dt(ChVector<>(0,0,0));
	shipPtr->SetMass(boxMass);
	shipPtr->SetInertiaXX(ChVector<>(bI2, bI3, bI1));
	shipPtr->SetIdentifier(shipId);
	shipPtr->SetCollide(true);
	shipPtr->SetBodyFixed(false);

	shipPtr->GetCollisionModel()->ClearModel();
//	shipPtr->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelop); //envelop is .03 by default
	utils::AddBoxGeometry(shipPtr.get_ptr(), 0.5 * ChVector<>(box_X, box_Y, box_Z), ChVector<>(0,0,0)); //beginning wall 2. Need "0.5 *" since chronoparallel is apparently different
	shipPtr->GetCollisionModel()->BuildModel();
	mphysicalSystem.Add(shipPtr);
}

//***** prismatic constraint between ship and bed
void Add_ship_ground_prismatic(ChSystemParallelDVI& mphysicalSystem) {
	ChSharedPtr<ChLinkLockPrismatic> shipGroundPrismatic(new ChLinkLockPrismatic);
	shipGroundPrismatic->Initialize(shipPtr, bin,
			ChCoordsys<>(ChVector<>(.30,  .09, -.25) , QUNIT)
			);
	shipGroundPrismatic->SetName("ship_ground_prismatic");
	mphysicalSystem.AddLink(shipGroundPrismatic);
}

void MoveShip_Kinematic(double dT) {
	ChVector<> shipVel = ChVector<>(0,0,shipVelocity);
	ChVector<> shipPos = shipInitialPos + shipVel * (dT - timePause);
	shipPtr->SetPos(shipPos);
	shipPtr->SetPos_dt(ChVector<>(0,0,shipVelocity));
	shipPtr->SetRot(ChQuaternion<>(1,0,0,0));
	shipPtr->SetWvel_loc(ChVector<>(0,0,0));
}

void FixShip_Kinematic() {
	shipPtr->SetPos_dt(ChVector<>(0,0,0));
}

void Add_Actuator(ChSystemParallelDVI& mphysicalSystem) {
//	shipPtr->SetPos_dt(ChVector<>(0,0,shipVelocity));
	ChSharedPtr<ChLinkLinActuator> actuator(new ChLinkLinActuator);
	ChVector<> pt1 = shipPtr->GetPos();
	ChVector<> pt2 = pt1 + ChVector<>(0, 0, 1000); //a large number in the z direction
	actuator->Initialize(shipPtr, bin, false, ChCoordsys<>(pt1, QUNIT), ChCoordsys<>(pt2, QUNIT));
	actuator->SetName("actuator");
	actuator->Set_lin_offset((pt2 - pt1).Length());
	mphysicalSystem.AddLink(actuator);
	//*** avoid infinite acceleration


}

void MoveShip_Actuator(
		ChSystemParallelDVI& mphysicalSystem,
		ChSharedPtr<ChFunction> actuator_fun,
		int path_piece) {
	static int interval = -1;
	if (path_piece == interval) return;
	ChSharedPtr<ChLinkLinActuator> actuator;
	actuator = mphysicalSystem.SearchLink("actuator").StaticCastTo<ChLinkLinActuator>();
	actuator->Set_dist_funct(actuator_fun);
	interval = path_piece;
}

//void MoveShip_PID(ChSystemParallelDVI& mphysicalSystem) {
//	char forceTag[] = "pulling_force";
//	ChSharedPtr<ChControllerPID> my_controllerPID;
//	static bool onCall = false;
//	if (!onCall) {
//		onCall = true;
//		my_controllerPID = ChSharedPtr<ChControllerPID>(new ChControllerPID);
//		my_controllerPID->P = 1.0e9;
//		my_controllerPID->D = 1.0e8;
//		my_controllerPID->I = 1.0e8;
//		ChSharedPtr<ChForce> pullingForce = ChSharedPtr<ChForce>(new ChForce);
//		pullingForce->SetMode(FTYPE_FORCE); // no need for this. It is the default option.
//		shipPtr->AddForce(pullingForce);
//		// ** or: hydroForce = ChSharedPtr<ChForce>(new ChForce());
//		pullingForce->SetName(forceTag);
//		pullingForce->SetVpoint(shipPtr->GetPos());
//		pullingForce->SetMforce(0);
//		pullingForce->SetDir(ChVector<>(1,0,0));
//	}
//	// Arman: check PID syntax. Do you need some velocity stuff?
//	double forcePID_X = my_controllerPID->Get_Out(shipPtr->GetBody()->GetPos().z - shipInitialPosZ - shipVelocity * (mphysicalSystem.GetChTime() - timePause), mphysicalSystem.GetChTime());
//	ChSharedPtr<ChForce> pullingForce = shipPtr->GetBody()->SearchForce(forceTag);
//	pullingForce->SetMforce(forcePID_X);
//	pullingForce->SetDir(ChVector<>(0,0,-1));
//}
 
int main(int argc, char* argv[])
{ 
	ChTimer<double> myTimerTotal;
	ChTimer<double> myTimerStep;
	int threads = 2;
	MySeed(964);

	myTimerTotal.start();
	ofstream outSimulationInfo("SimInfo.txt");

	if (argc > 1) {
		const char* text = argv[1];
    	threads = atoi(text);
	}
	outSimulationInfo << "** num threads: " << threads << endl;

	// ***** params
	double gravity = 9.81;
	double dT = 0.1* mradius / shipVelocity; //moving 0.1*R at each time step
	double time_end = 100;
	double out_fps = 50;
	uint max_iteration = 5000;//10000;
	double tolerance = 1e-3;
	// ************


	// Create a ChronoENGINE physical system
	ChSystemParallelDVI mphysicalSystem;

	//******************** OMP settings **************
	// Set number of threads.
	int max_threads = mphysicalSystem.GetParallelThreadNumber();
	outSimulationInfo << "** max number of threads: " << max_threads << endl;
	if (threads > max_threads)
	  threads = max_threads;
	mphysicalSystem.SetParallelThreadNumber(threads);
	omp_set_num_threads(threads);
	//************************************************

	// Set gravitational acceleration
	mphysicalSystem.Set_G_acc(ChVector<>(0, -gravity, 0));

	// Set solver parameters
	mphysicalSystem.GetSettings()->solver.solver_mode = SLIDING; //NORMAL, SPINNING
	mphysicalSystem.GetSettings()->solver.max_iteration_normal = max_iteration / 3;
	mphysicalSystem.GetSettings()->solver.max_iteration_sliding = max_iteration / 3;
	mphysicalSystem.GetSettings()->solver.max_iteration_spinning = 0;
	mphysicalSystem.GetSettings()->solver.max_iteration_bilateral = max_iteration / 3;
	mphysicalSystem.GetSettings()->solver.tolerance = 0;//tolerance;
	mphysicalSystem.GetSettings()->solver.alpha = 0;  //Arman, find out what is this
	mphysicalSystem.GetSettings()->solver.contact_recovery_speed = 2 * shipVelocity;  //Arman, I hope it is the counterpart of SetMaxPenetrationRecoverySpeed
	mphysicalSystem.ChangeSolverType(APGDRS);  //Arman check this APGD APGDBLAZE
	mphysicalSystem.GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;

	mphysicalSystem.GetSettings()->collision.collision_envelope = collisionEnvelop;
	mphysicalSystem.GetSettings()->collision.bins_per_axis = I3(10, 10, 10); //Arman check
//	mphysicalSystem.GetSettings()->collision.min_body_per_bin = 50;		// Arman check: According to Hammad these are not used anymore
//	mphysicalSystem.GetSettings()->collision.max_body_per_bin = 100;	// Arman check: According to Hammad these are not used anymore

	//******************* Irrlicht and driver types **************************
#define irrlichtVisualization true
	driveType = ACTUATOR;//KINEMATIC : ACTUATOR
	//************************************************************************
	outSimulationInfo << "****************************************************************************" << endl;
	outSimulationInfo << "dT: " << dT <<" shipVelocity: "<< shipVelocity << " particles_radius: " << mradius <<
			" timePause: " << timePause << " timeMove: " << timeMove << endl;

	ofstream outForceData("forceData.txt");

	// Create all the rigid bodies.
	create_system_particles(mphysicalSystem);
	Add_ship_ground_prismatic(mphysicalSystem);

#ifdef CHRONO_PARALLEL_HAS_OPENGL2
   opengl::ChOpenGLWindow &gl_window = opengl::ChOpenGLWindow::getInstance();
   gl_window.Initialize(1280, 720, "mixerDVI", &mphysicalSystem);
   gl_window.SetCamera(ChVector<>(-3,12,-8), ChVector<>(7.2, 6, 8.2), ChVector<>(0, 1, 0)); //camera

   // Uncomment the following two lines for the OpenGL manager to automatically
   // run the simulation in an infinite loop.
   //gl_window.StartDrawLoop(time_step);
   //return 0;
#endif

#if irrlichtVisualization
		// Create the Irrlicht visualization (open the Irrlicht device,
		// bind a simple user interface, etc. etc.)
		ChIrrApp application(&mphysicalSystem, L"Bricks test",core::dimension2d<u32>(800,600),false, true);
		// Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
		ChIrrWizard::add_typical_Logo  (application.GetDevice());
		ChIrrWizard::add_typical_Sky   (application.GetDevice());
		ChIrrWizard::add_typical_Lights(application.GetDevice(), core::vector3df(14.0f, 44.0f, -18.0f), core::vector3df(-3.0f, 8.0f, 6.0f), 59,  40);
		ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(7.2,20,21.2), core::vector3df(7.2,1,10)); //   (7.2,30,0) :  (-3,12,-8)
		// Use this function for adding a ChIrrNodeAsset to all items
		// If you need a finer control on which item really needs a visualization proxy in
		// Irrlicht, just use application.AssetBind(myitem); on a per-item basis.
		application.AssetBindAll();
		// Use this function for 'converting' into Irrlicht meshes the assets
		// into Irrlicht-visualizable meshes
		application.AssetUpdateAll();

		application.SetStepManage(true);
		application.SetTimestep(dT);  					//Arman modify
#endif
		double iceThickness = numLayers * 2 * mradius * cos(CH_C_PI / 6.0);
		double calculatedPorosity = 1 - (mphysicalSystem.Get_bodylist()->size() - 2) * 4.0 / 3 *CH_C_PI * pow(mradius, 3) /
				( hdim.x * hdim.z * iceThickness);
		outForceData << "[1] time, [2-5] forceContact (x, y, z, magnitude), [6-9] forceActuator (x, y, z, magnitude), [10-13] Ice pressure contact (x, y, z, magnitude), [14-17] Ice pressure actuator (x, y, z, magnitude), [18] shipPos, [19] shipVel, [20] energy, [21] iceThickness, [22] timePerStep, [23] timeElapsed. ## numSpheres" << mphysicalSystem.Get_bodylist()->end() - mphysicalSystem.Get_bodylist()->begin()
				<< " pauseTime: " << timePause<< " setVelocity: "<< shipVelocity << " ship_width: " << ship_width << " set porosity: " << porosity << " calculated porosity: " << calculatedPorosity << endl;
		outForceData.close();
		outSimulationInfo << "Real Time, Compute Time" << endl;

	outSimulationInfo << "***** number of bodies: " << mphysicalSystem.Get_bodylist()->size() << endl;
	bool moveTime = false;
	//****************************************** Time Loop *************************************
	ChSharedPtr<ChFunction_Const> actuator_fun0(new ChFunction_Const(0));
	ChSharedPtr<ChFunction_Ramp> actuator_fun1(new ChFunction_Ramp(0, -shipVelocity));
	if (driveType == ACTUATOR) {
		Add_Actuator(mphysicalSystem);
	}

	while(mphysicalSystem.GetChTime() < timeMove+timePause) //arman modify
	{
		myTimerStep.start();
		// ****** include force or motion ********
		switch (driveType) {
		case KINEMATIC:
			(mphysicalSystem.GetChTime() < timePause) ?
				FixShip_Kinematic() :
				MoveShip_Kinematic(mphysicalSystem.GetChTime());
			break;
		case ACTUATOR:
			(mphysicalSystem.GetChTime() < timePause) ?
				MoveShip_Actuator(mphysicalSystem, actuator_fun0.StaticCastTo<ChFunction>(), 0) :
				MoveShip_Actuator(mphysicalSystem, actuator_fun1.StaticCastTo<ChFunction>(), 1);
		}
		// ****** end of force or motion *********
#if irrlichtVisualization
		if ( !(application.GetDevice()->run()) ) break;
		application.GetVideoDriver()->beginScene(true, true, SColor(255,140,161,192));
		ChIrrTools::drawGrid(application.GetVideoDriver(), 1,1, 40,40,
			ChCoordsys<>(ChVector<>(0,-2,0),Q_from_AngAxis(CH_C_PI/2,VECT_X)), video::SColor(50,90,90,150),true);
		application.DrawAll();
		application.DoStep();
		application.GetVideoDriver()->endScene();
#else
#ifdef CHRONO_PARALLEL_HAS_OPENGL2
		if (gl_window.Active()) {
		 gl_window.DoStepDynamics(dT);
		 gl_window.Render();
		}
#else
		mphysicalSystem.DoStepDynamics(dT);
#endif
#endif

		//******************** ship force*********************
		ChVector<> mForceActuator = ChVector<>(0,0,0);
		ChVector<> mTorqueActuator = ChVector<>(0,0,0);
		ChVector<> icePressureActuator = ChVector<>(0,0,0);
		ChVector<> mTorqueContact;
		ChVector<> mForceContact;


		if (driveType == ACTUATOR) {
			ChSharedPtr<ChLinkLinActuator> actuator;
			actuator = mphysicalSystem.SearchLink("actuator").StaticCastTo<ChLinkLinActuator>();
			mForceActuator = actuator->Get_react_force();
			mTorqueActuator = actuator->Get_react_torque();
			icePressureActuator = mForceActuator / iceThickness / ship_width;
		}
		calc_ship_contact_forces(mphysicalSystem, mForceContact, mTorqueContact);
		ChVector<> icePressureContact = mForceContact / iceThickness / ship_width;

		myTimerStep.stop();
		myTimerTotal.stop();
		//****************************************************
		vector<ChBody*>::iterator ibody = mphysicalSystem.Get_bodylist()->begin();
		double energy = 0;
		while (ibody != mphysicalSystem.Get_bodylist()->end()) {
			create_hydronynamic_force(*ibody, mphysicalSystem, surfaceLoc, false);
			energy += pow((*ibody)->GetPos_dt().Length() , 2);
			ibody++;
		}

		printf("*** total number of contacts %d, num bodies %d\n", mphysicalSystem.GetNcontacts(), mphysicalSystem.Get_bodylist()->size());
		stringstream outDataSS;
		outDataSS << mphysicalSystem.GetChTime() << ", " <<
				mForceContact.x << ", " << mForceContact.y << ", " << mForceContact.z << ", " << mForceContact.Length() << ", " <<
				mForceActuator.x << ", " << mForceActuator.y << ", " << mForceActuator.z << ", " << mForceActuator.Length() << ", " <<
				icePressureContact.x << ", " << icePressureContact.y << ", " << icePressureContact.z << ", " << icePressureContact.Length() << ", " <<
				icePressureActuator.x << ", " << icePressureActuator.y << ", " << icePressureActuator.z << ", " << icePressureActuator.Length() << ", " <<
				shipPtr->GetPos().z << ", " << shipPtr->GetPos_dt().z << ", " << energy << ", " << iceThickness  << ", " << myTimerStep() << ", " << myTimerTotal() << endl;
		ofstream outData("forceData.txt", ios::app);
		outData<<outDataSS.str();
		outData.close();

		double numIter = ((ChLcpSolverParallelDVI*)mphysicalSystem.GetLcpSolverSpeed())->GetTotalIterations();
		outSimulationInfo << "Time: " <<  mphysicalSystem.GetChTime() <<
				" Ship pos: " << shipPtr->GetPos().x << ", " << shipPtr->GetPos().y << ", " <<  shipPtr->GetPos().z <<
				" Ship vel: " << shipPtr->GetPos_dt().x << ", " << shipPtr->GetPos_dt().y << ", " <<  shipPtr->GetPos_dt().z <<
				" energy: " << energy <<
				" time per step: " << myTimerStep() <<
				" time elapsed: " << myTimerTotal() <<
				" Ship force: " << mForceContact.x << ", " << mForceContact.y << ", " <<  mForceContact.z <<
				" ice thickness: " << iceThickness <<
				" number of Iteration: " << numIter << endl;
		cout << "Time: " <<  mphysicalSystem.GetChTime() <<
				" Ship vel: " << shipPtr->GetPos_dt().x << ", " << shipPtr->GetPos_dt().y << ", " <<  shipPtr->GetPos_dt().z <<
				" energy: " << energy <<
				" time per step: " << myTimerStep() <<
				" time elapsed: " << myTimerTotal() <<
				" Ship force: " << mForceContact.x << ", " << mForceContact.y << ", " <<  mForceContact.z <<
				" ice thickness: " << iceThickness <<
				" number of Iteration: " << numIter << endl;

	}

	outForceData.close();
	return 0;
}
  
