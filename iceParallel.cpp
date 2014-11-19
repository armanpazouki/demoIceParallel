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
///////////////////////////////////////////////////
    
  
 
////#include "physics/ChBodyEasy.h"
//#include "physics/ChContactContainer.h"
//#include "collision/ChCModelBulletBody.h"
////#include "core/ChTimer.h"
////#include "core/ChRealtimeStep.h"
////#include "assets/ChTexture.h"
////#include "unit_IRRLICHT/ChIrrApp.h"
//#include <cstring>
//#include <fstream>
////#include <map>

//*************** chrono parallel
#include <stdio.h>
#include <vector>
#include <cmath>

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

#include "chrono_utils/ChUtilsCreators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

#ifdef CHRONO_PARALLEL_HAS_OPENGL2
#include "chrono_opengl/ChOpenGLWindow.h"
#endif
//***********************************
// Use the namespace of Chrono

using namespace chrono;
using namespace chrono::collision;

// Use the main namespaces of Irrlicht
//using namespace irr;
//using namespace core;
//using namespace scene;
//using namespace video;
//using namespace io;
//using namespace gui;
//using namespace std;

const double rhoF = 1000;
const double rhoR = 917;
const double mu_Viscosity = .001;//.1;
const ChVector<> surfaceLoc = ChVector<>(0, .04, -.08);

//******************* ship and sphere stuff
double mradius = .4;
int numLayers = 2;

//ChBodySceneNode* shipPtr;
//ChSharedPtr<ChBodyEasyBox> shipPtr; //old irrlicht version
ChSharedBodyPtr shipPtr;
const double shipVelocity = 5.4;//.27;//1; //arman modify
double shipInitialPosZ = 0;
const double timePause = 1.0; //arman modify
const double timeMove = 2.5;
double ship_width = 4;
double box_X = ship_width, box_Y = 10, box_Z = .4;
double collisionEnvelop = .04 * mradius;
ChVector<> shipInitialPos;
//**********************************

void Calc_Hydrodynamics_Forces(ChVector<> & F_Hydro, ChVector<> & forceLoc, ChVector<> & T_Drag,
		ChBody* mrigidBody, ChSystemParallelDVI& mphysicalSystem, const chrono::ChVector<>& freeSurfaceLocation) {
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

void create_hydronynamic_force(ChBody* mrigidBody, ChSystemParallelDVI& mphysicalSystem, const chrono::ChVector<>& freeSurfaceLocation, bool createForce) {
	// ***** insertion of force
	ChSharedPtr<ChForce> hydroForce;
	ChSharedPtr<ChForce> hydroTorque;

//	std::string forceTag("hydrodynamics_force");
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

void calc_ship_contact_forces(ChSystemParallelDVI& mphysicalSystem, ChVector<> & mForce, ChVector<> & mTorque) {
	mForce = ChVector<>(0,0,0);
	mTorque = ChVector<>(0,0,0);
	ChContactContainer* container  = (ChContactContainer *) mphysicalSystem.GetContactContainer();
//	std::map<ChBody*, ChVector<> > m_forces;
//	std::map<ChBody*, ChVector<> > m_torques;
//	ChVector<> mForce;
//	ChVector<> mTorque;

	std::list<ChContact*> m_list = container->GetContactList();
	for (std::list<ChContact *>::iterator it=m_list.begin(); it != m_list.end(); ++it){
	  ChVector<> force_contactFrame = (*it)->GetContactForce();
	  ChVector<> force_abs = *((*it)->GetContactPlane()) * force_contactFrame;
	  ChModelBulletBody * model_A = (ChModelBulletBody *) (*it)->GetModelA();
	  ChModelBulletBody * model_B = (ChModelBulletBody *) (*it)->GetModelB();
//		  if ((model_A->GetBody() != shipPtr->GetBody()) && (model_B->GetBody() != shipPtr->GetBody())) {
//			  continue;
//		  }
	  ChBody * body_A = model_A->GetBody();
	  ChBody * body_B = model_B->GetBody();

	  if (body_A == (ChBody*)shipPtr.get_ptr()) {
		  mForce -= force_abs;

		  ChVector<> point_on_A = (*it)->GetContactP1();
		  mTorque -= (point_on_A - body_A->GetPos()) % force_abs;
//		  ChVector<> local_point_on_A = ChTransform<>::TransformParentToLocal(point_on_A, body_A->GetPos(), body_A->GetRot());
//		  mTorque += local_point_on_A % force;
	  } else if (body_B == (ChBody*)shipPtr.get_ptr()) {
		  mForce += force_abs;

		  ChVector<> point_on_B = (*it)->GetContactP2();
		  mTorque += (point_on_B - body_B->GetPos()) % force_abs;
//		  ChVector<> local_point_on_B = ChTransform<>::TransformParentToLocal(point_on_B, body_B->GetPos(), body_B->GetRot());
//		  mTorque -= local_point_on_B % force;
	  }
	}
}
//***********************************
void CreateSphere(ChSystemParallelDVI& mphysicalSystem, ChSharedBodyPtr mrigidBody, ChVector<> pos, double mmass, double minert) {
	static int sphereID = 0;

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
//		mrigidBody->GetMaterialSurface()->SetFriction(0.4f);
//		mrigidBody->GetMaterialSurface()->SetCompliance(0.0);
//		mrigidBody->GetMaterialSurface()->SetComplianceT(0.0);
//		mrigidBody->GetMaterialSurface()->SetDampingF(0.2);
//
//		mrigidBody->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelop); //envelop is .03 by default
	mrigidBody->GetCollisionModel()->ClearModel();
	utils::AddSphereGeometry(mrigidBody.get_ptr(), mradius);
	mrigidBody->GetCollisionModel()->BuildModel();

	create_hydronynamic_force(mrigidBody.get_ptr(), mphysicalSystem, surfaceLoc, true);  //Arman : hydrodynamic forces

//	// optional, attach a texture for better visualization  //Arman irrlicht
//	ChSharedPtr<ChTexture> mtextureball(new ChTexture());
//	mtextureball->SetTextureFilename(GetChronoDataFile("../data/bluwhite.png"));
//	mrigidBody->AddAsset(mtextureball);

	mphysicalSystem.Add(mrigidBody);
//    sys->AddBody(mrigidBody);
}
//***********************************
void GenerateIceLayers_Rectangular(
		ChSystemParallelDVI& mphysicalSystem,
		ChVector<> boxMin,
		ChVector<> boxMax,
		double global_x, //global offset in x
		double global_y, //global offset in y
		double global_z,
		double expandR,
		double mmass,
		double minert) //global offset in z
{
	int numColX = (boxMax.x - boxMin.x - expandR) / (2 * expandR);
	int numColZ = (boxMax.z - boxMin.z - expandR) / (2 * expandR);
	double spacing = 2 * expandR;
	for (int j = 0; j < numLayers; j++) {
		for (int i = 0; i < numColX; i++) {
			for (int k = 0; k < numColZ; k++) {
				// Create a ball that will collide with wall
				ChSharedBodyPtr mrigidBody(new ChBody(new ChCollisionModelParallel));
				ChVector<> pos = ChVector<>(
						ChVector<>(global_x, global_y, global_z)
						+ ChVector<>((i+0.5) * spacing, j * spacing, (k+0.5) * spacing)
						+ .5 * (spacing - 2 * mradius) * ChVector<>(ChRandom(), ChRandom(), ChRandom())
						);
				CreateSphere(mphysicalSystem, mrigidBody, pos, mmass, minert);
			}
		}

	}
}
//***********************************
void addHCPSheet(
		ChSystemParallelDVI& mphysicalSystem,
		int grid_x,       //number of particles in x direction
		int grid_z,       //number of particles in z direction
		double height,    //height of layer
		double global_x,  //global offset of sheet in x
		double global_z,
		double expandR,
		double mmass,
		double minert)  //global offset of sheet in z
{
    double offset = 0;
    double x = 0, y = height, z = 0;
    for (int i = 0; i < grid_x; i++) {
      for (int k = 0; k < grid_z; k++) {
        //need to offset alternate rows by radius
        offset = (k % 2 != 0) ? expandR : 0;
        //x position, shifted to center
        x = i * 2 * expandR + offset + expandR + global_x;
        //z position shifted to center
        z = k * (sqrt(3.0) * expandR)  + expandR + global_z;
        // x, y, z contain coordinates for sphere position





		ChSharedBodyPtr mrigidBody(new ChBody(new ChCollisionModelParallel));
		ChVector<> pos = ChVector<>(x, height, z);
		CreateSphere(mphysicalSystem, mrigidBody, pos, mmass, minert);
      }
    }
}
//***********************************
void GenerateIceLayers_Hexagonal(
		ChSystemParallelDVI& mphysicalSystem,
		ChVector<> boxMin,
		ChVector<> boxMax,
		double global_x, //global offset in x
		double global_y, //global offset in y
		double global_z,
		double expandR,
		double mmass,
		double minert) //global offset in z
{
	int numColX = (boxMax.x - boxMin.x - expandR) / (2 * expandR);
	int numColZ = (boxMax.z - boxMin.z - expandR) / (sqrt(3.0) * expandR);
    double offset_x = 0, offset_z = 0, height = 0;
    for (int j = 0; j < numLayers; j++) {
      height = j * (sqrt(3.0) * expandR);
      //need to offset each alternate layer by radius in both x and z direction
      offset_x = offset_z = (j % 2 != 0) ? expandR : 0;
      addHCPSheet(mphysicalSystem, numColX, numColZ, height + global_y, offset_x+global_x, offset_z+global_z, expandR, mmass, minert);
    }
}
//***********************************
void create_ice_particles(ChSystemParallelDVI& mphysicalSystem)
{
	//**************** sphere prob
	double expandR = mradius*1.05;

	double iceThickness = numLayers * mradius * 2;
	double buttomLayerDY = rhoR / rhoF *  iceThickness - mradius;

	double mmass = (4./3.)*CH_C_PI*pow(mradius,3)*rhoR;
	double minert = (2./5.)* mmass * pow(mradius,2);
	ChVector<> boxMin = ChVector<>(-.8, 0, -2.4);
	ChVector<> hdim = ChVector<>(16, 12, 21.2);
	ChVector<> centerLF = 0.5 * hdim;
	ChVector<> center = 0.5 * hdim + boxMin;
	ChVector<> boxMax = boxMin + hdim;
	printf("************************** Generate Ice, ButtomLayer_Y %f\n", buttomLayerDY);
	double global_x = boxMin.x;
	double global_y = surfaceLoc.y - buttomLayerDY;
	double global_z = boxMin.z;
//	GenerateIceLayers_Rectangular(mphysicalSystem,
//			boxMin, boxMax,
//			global_x, global_y, global_z,
//			expandR, mmass,	minert);
	GenerateIceLayers_Hexagonal(mphysicalSystem,
			boxMin, boxMax,
			global_x, global_y, global_z,
			expandR, mmass,	minert);

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

	ChSharedBodyPtr bin(new ChBody(new ChCollisionModelParallel));
	bin->SetMaterialSurface(mat);
	bin->SetIdentifier(binId);
	bin->SetMass(1);
	bin->SetPos(center);
	bin->SetRot(ChQuaternion<>(1, 0, 0, 0));
	bin->SetCollide(true);
	bin->SetBodyFixed(true);

	bin->GetCollisionModel()->ClearModel();
	//utils::AddBoxGeometry(bin.get_ptr(), ChVector<>(110,1,110), ChVector<>(0, -10, 0)); //earth, not necessary

	utils::AddBoxGeometry(bin.get_ptr(), ChVector<>(hdim.x, hdim.y, hthick), boxMin + ChVector<>(centerLF.x, 0, hdim.z + 0.5*hthick));	//end wall
	utils::AddBoxGeometry(bin.get_ptr(), ChVector<>(hthick, hdim.y, hdim.z), boxMin + ChVector<>(0 - 0.5*hthick, 0, centerLF.z));		//side wall
	utils::AddBoxGeometry(bin.get_ptr(), ChVector<>(hthick, hdim.y, hdim.z), boxMin + ChVector<>(hdim.x + 0.5*hthick, 0, centerLF.z));	//side wall
	utils::AddBoxGeometry(bin.get_ptr(), ChVector<>(small_wall_Length, hdim.y, hthick), boxMin + ChVector<>(0.5*small_wall_Length, 0, -0.5*hthick)); 	//beginning wall 1
	utils::AddBoxGeometry(bin.get_ptr(), ChVector<>(small_wall_Length, hdim.y, hthick), boxMin + ChVector<>(hdim.x - 0.5*small_wall_Length, 0, -0.5*hthick)); //beginning wall 2
	bin->GetCollisionModel()->BuildModel();

	mphysicalSystem.AddBody(bin);

	//**************** create ship
	double boxMass = rhoR * box_X * box_Y * box_Z;
	double bI1 = 1.0 / 12 * boxMass * (pow(box_X, 2) + pow(box_Y, 2));
	double bI2 = 1.0 / 12 * boxMass * (pow(box_Y, 2) + pow(box_Z, 2));
	double bI3 = 1.0 / 12 * boxMass * (pow(box_X, 2) + pow(box_Z, 2));
	shipInitialPosZ = boxMin.z - .5 * box_Z;

	shipPtr = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	shipInitialPos = ChVector<>(center.x,  1, shipInitialPosZ);
	shipPtr->SetPos(shipInitialPos);
	shipPtr->SetRot(ChQuaternion<>(1,0,0,0));
	shipPtr->SetMaterialSurface(mat);
	shipPtr->SetPos_dt(ChVector<>(0,0,0));
	shipPtr->SetMass(boxMass);
	shipPtr->SetInertiaXX(ChVector<>(bI2, bI3, bI1));
	bin->SetIdentifier(binId);
	bin->SetCollide(true);
	bin->SetBodyFixed(false);

	bin->GetCollisionModel()->ClearModel();
//	shipPtr->GetCollisionModel()->SetDefaultSuggestedEnvelope(collisionEnvelop); //envelop is .03 by default
	utils::AddBoxGeometry(bin.get_ptr(), ChVector<>(box_X, box_Y, box_Z), shipInitialPos); //beginning wall 2
	bin->GetCollisionModel()->BuildModel();
	mphysicalSystem.Add(shipPtr);

//	// optional, attach a texture for better visualization, Arman Texture
//	ChSharedPtr<ChTexture> mtexturebox(new ChTexture());
//	mtexturebox->SetTextureFilename(GetChronoDataFile("../data/cubetexture_borders.png"));
//	shipPtr->AddAsset(mtexturebox);
//	//**** end ship initialization

	char forceTag[] = "pulling_force";
	ChSharedPtr<ChForce> pullingForce = ChSharedPtr<ChForce>(new ChForce);
	pullingForce->SetMode(FTYPE_FORCE); // no need for this. It is the default option.
	shipPtr->AddForce(pullingForce);
	// ** or: hydroForce = ChSharedPtr<ChForce>(new ChForce());
	pullingForce->SetName(forceTag);
	pullingForce->SetVpoint(shipPtr->GetPos());
	pullingForce->SetMforce(0);
	pullingForce->SetDir(ChVector<>(1,0,0));

	//***** prismatic constraint between ship and bed
//	ChSharedPtr<ChLinkLockPlanePlane> shipConstraint(new ChLinkLockPlanePlane);
//	shipConstraint->Initialize(shipPtr->GetBody(), earthPtr->GetBody(),
//			ChCoordsys<>(ChVector<>(30,  9, -25) , Q_from_AngAxis(CH_C_PI/2, VECT_X))
//			);
	ChSharedPtr<ChLinkLockPrismatic> shipConstraint(new ChLinkLockPrismatic);
	shipConstraint->Initialize(shipPtr, bin,
			ChCoordsys<>(ChVector<>(.30,  .09, -.25) , QUNIT)
			);
	mphysicalSystem.AddLink(shipConstraint);
}

void MoveShip(ChSystemParallelDVI& mphysicalSystem) {
	static bool onCall = false;
//	if (!onCall) {
//		onCall = true;
		ChVector<> shipVel = ChVector<>(0,0,shipVelocity);
		ChVector<> shipPos = shipInitialPos + shipVel * (mphysicalSystem.GetChTime() - timePause);
		shipPtr->SetPos(shipPos);
		shipPtr->SetPos_dt(ChVector<>(0,0,shipVelocity));
		shipPtr->SetRot(ChQuaternion<>(1,0,0,0));
		shipPtr->SetWvel_loc(ChVector<>(0,0,0));


//	}
//    ChSharedPtr<ChControllerPID> my_controllerPID(new ChControllerPID);
//    my_controllerPID->P = 1.0e9;
//    my_controllerPID->D = 1.0e8;
//    my_controllerPID->I = 1.0e8;
//
//    double forcePID_X = my_controllerPID->Get_Out(shipPtr->GetBody()->GetPos().z - shipInitialPosZ - shipVelocity * (mphysicalSystem.GetChTime() - timePause), mphysicalSystem.GetChTime());
//    char forceTag[] = "pulling_force";
//	ChSharedPtr<ChForce> pullingForce = shipPtr->GetBody()->SearchForce(forceTag);
//	pullingForce->SetMforce(forcePID_X);
//	pullingForce->SetDir(ChVector<>(0,0,-1));
}

void FixShip(ChSystemParallelDVI& mphysicalSystem) {
	shipPtr->SetPos_dt(ChVector<>(0,0,0));
}
 
int main(int argc, char* argv[])
{ 
	ChTimer<double> myTimer;
	int threads = 8;

	// ***** params
	double gravity = 9.81;
	double dT = 0.1* mradius / shipVelocity; //moving 0.1*R at each time step
//	double time_step = 1e-3;
	double time_end = 100;
	double out_fps = 50;
	uint max_iteration = 30;
	real tolerance = 1e-3;
	// ************

#define irrlichtVisualization false
	// Create a ChronoENGINE physical system
	ChSystemParallelDVI mphysicalSystem;

	//******************** OMP settings **************
	// Set number of threads.
	int max_threads = mphysicalSystem.GetParallelThreadNumber();
	if (threads > max_threads)
	  threads = max_threads;
	mphysicalSystem.SetParallelThreadNumber(threads);
	omp_set_num_threads(threads);
	//************************************************
	// Set gravitational acceleration
	mphysicalSystem.Set_G_acc(ChVector<>(0, -gravity, 0));

	// Set solver parameters
	mphysicalSystem.GetSettings()->solver.solver_mode = SLIDING; //Arman, What is this
	mphysicalSystem.GetSettings()->solver.max_iteration_normal = max_iteration / 3;
	mphysicalSystem.GetSettings()->solver.max_iteration_sliding = max_iteration / 3;
	mphysicalSystem.GetSettings()->solver.max_iteration_spinning = 0;
	mphysicalSystem.GetSettings()->solver.max_iteration_bilateral = max_iteration / 3;
	mphysicalSystem.GetSettings()->solver.tolerance = tolerance;
	mphysicalSystem.GetSettings()->solver.alpha = 0;  //Arman, find out what is this
	mphysicalSystem.GetSettings()->solver.contact_recovery_speed = 2 * shipVelocity;  //Arman, I hope it is the counterpart of SetMaxPenetrationRecoverySpeed
	mphysicalSystem.ChangeSolverType(APGDRS);  //Arman check this
	mphysicalSystem.GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;

	mphysicalSystem.GetSettings()->collision.collision_envelope = collisionEnvelop;
	mphysicalSystem.GetSettings()->collision.bins_per_axis = I3(10, 10, 10); //Arman check
	mphysicalSystem.GetSettings()->collision.min_body_per_bin = 50;			// Arman check
	mphysicalSystem.GetSettings()->collision.max_body_per_bin = 100;		// Arman check

	printf("****************************************************************************\n");
	printf("dT: %f, shipVelocity: %f, particles_radius: %f, timePause: %f, timeMove: %f\n\n", dT, shipVelocity, mradius, timePause, timeMove);



	std::fstream outForceData("forceData.txt", std::ios::out);

	// Create all the rigid bodies.
	create_ice_particles(mphysicalSystem);

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
		ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(-3,12,-8), core::vector3df(0,1,0));
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

	outForceData << "time, forceX, forceY, forceZ, forceMag, pressureX, pressureY, pressureZ, pressureMag, shipVelocity, shipPosition, energy, timePerStep.## numSpheres" << mphysicalSystem.Get_bodylist()->end() - mphysicalSystem.Get_bodylist()->begin()
			<< " pauseTime: " << timePause<< " setVelocity: "<< shipVelocity << std::endl;

	while(mphysicalSystem.GetChTime() < timeMove+timePause) //arman modify
	{
		myTimer.start();
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
		if (mphysicalSystem.GetChTime() > timePause) {
			MoveShip(mphysicalSystem);
		} else {
			FixShip(mphysicalSystem);
		}
		//******************** ship force*********************
//		ChVector<> shipForce = shipPtr->GetBody()->Get_Xforce();
//		printf("force %f\n",shipForce.z);

		ChVector<> mForce;
		ChVector<> mTorque;
		calc_ship_contact_forces(mphysicalSystem, mForce, mTorque);
		ChVector<> icePressure = mForce / (numLayers * 2 * mradius * cos(CH_C_PI / 6)) / ship_width;


		myTimer.stop();
		//****************************************************
//		for(int i=0; i<mphysicalSystem.Get_bodylist()->size(); i++){
//			create_hydronynamic_force(mphysicalSystem.Get_bodylist()->at(i), mphysicalSystem, surfaceLoc, false);
//
//		}
		std::vector<ChBody*>::iterator ibody = mphysicalSystem.Get_bodylist()->begin();
		double energy = 0;
		while (ibody != mphysicalSystem.Get_bodylist()->end()) {
			create_hydronynamic_force(*ibody, mphysicalSystem, surfaceLoc, false);
			energy += pow((*ibody)->GetPos_dt().Length() , 2);
			ibody++;
		}
//		printf("time %f, force %f %f %f, shipVelocity %f, simulation time %f, energy %f\n", mphysicalSystem.GetChTime(), mForce.x, mForce.y, mForce.z, shipPtr->GetBody()->GetPos_dt().z, myTimer(), energy);
		outForceData << mphysicalSystem.GetChTime() << ", " << mForce.x << ", " << mForce.y << ", " << mForce.z << ", " <<
				mForce.Length() << ", " <<
				icePressure.x << ", " << icePressure.y << ", " << icePressure.z << ", " << icePressure.Length() << ", " <<
				shipPtr->GetPos_dt().z << ", " << shipPtr->GetPos().z << ", " << energy << ", " << myTimer() << std::endl;

		printf("Time %f, energy %f, time per step %f\n", mphysicalSystem.GetChTime(), energy, myTimer());
	}
	outForceData.close();
	return 0;
}
  
