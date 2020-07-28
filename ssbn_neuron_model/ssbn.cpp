/*
 *  ssbn.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include "ssbn.h"

#ifdef HAVE_GSL

// C++ includes:
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
//#include <random>
#include <sys/types.h>
#include <sys/stat.h>
//#include <stdlib.h>
#include <unistd.h>
// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include <fcntl.h>
/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< nest::ssbn >
  nest::ssbn::recordablesMap_;

namespace nest // template specialization must be placed in namespace
{
/*
 * Override the create() method with one call to RecordablesMap::insert_()
 * for each quantity to be recorded.
 */
template <>
void
RecordablesMap< ssbn >::create()
{
  // use standard names whereever you can for consistency!
  insert_(
    names::V_m, &ssbn::get_y_elem_< psdb::State_::V_M > );
  insert_( names::g_ex,
    &ssbn::get_y_elem_< psdb::State_::G_EXC > );
  insert_( names::g_in,
    &ssbn::get_y_elem_< psdb::State_::G_INH > );

  insert_( names::t_ref_remaining, &ssbn::get_r_ );
}
}

/* ----------------------------------------------------------------
 * Iteration function
 * ---------------------------------------------------------------- */

extern "C" inline int
nest::ssbn_dynamics( double,
  const double y[],
  double f[],
  void* pnode )
{
  // a shorthand
  typedef nest::ssbn::State_ S;

  // get access to node so we can almost work as in a member function
  assert( pnode );
  const nest::ssbn& node =
    *( reinterpret_cast< nest::ssbn* >( pnode ) );

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[].

  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away ...
  const double I_syn_exc = y[ S::G_EXC ] * ( y[ S::V_M ] - node.P_.E_ex );
  const double I_syn_inh = y[ S::G_INH ] * ( y[ S::V_M ] - node.P_.E_in );
  const double I_leak = node.P_.g_L * ( y[ S::V_M ] - node.P_.E_L );

  // dV_m/dt
  f[ 0 ] = ( -I_leak - I_syn_exc - I_syn_inh + node.B_.I_stim_ + node.P_.I_e )
    / node.P_.C_m;

  // d dg_exc/dt, dg_exc/dt
  f[ 1 ] = -y[ S::DG_EXC ] / node.P_.tau_synE;
  f[ 2 ] = y[ S::DG_EXC ] - ( y[ S::G_EXC ] / node.P_.tau_synE );

  // d dg_exc/dt, dg_exc/dt
  f[ 3 ] = -y[ S::DG_INH ] / node.P_.tau_synI;
  f[ 4 ] = y[ S::DG_INH ] - ( y[ S::G_INH ] / node.P_.tau_synI );

  return GSL_SUCCESS;
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::ssbn::Parameters_::Parameters_()
  : V_th( -55.0 )    // mV
  , V_reset( -60.0 ) // mV
  , t_ref( 2.0 )     // ms
  , g_L( 16.6667 )   // nS
  , C_m( 250.0 )     // pF
  , E_ex( 0.0 )      // mV
  , E_in( -85.0 )    // mV
  , E_L( -70.0 )     // mV
  , tau_synE( 0.2 )  // ms
  , tau_synI( 2.0 )  // ms
  , I_e( 0.0 )       // pA
  , n_s( 3 )  
{
}

nest::ssbn::State_::State_( const Parameters_& p )
  : r( 0 )
  , bm(0)
{
  y[ V_M ] = p.E_L; // initialize to reversal potential
  for ( size_t i = 1; i < STATE_VEC_SIZE; ++i )
    y[ i ] = 0;
}

nest::ssbn::State_::State_( const State_& s )
  : r( s.r ),
  bm(s.bm)
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
    y[ i ] = s.y[ i ];
}

nest::ssbn::State_& nest::psdb::State_::operator=(
  const State_& s )
{
  if ( this == &s ) // avoid assignment to self
    return *this;

  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
    y[ i ] = s.y[ i ];

  r = s.r;
  bm =  s.bm;
  return *this;
}

nest::ssbn::Buffers_::Buffers_( psdb& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

nest::ssbn::Buffers_::Buffers_( const Buffers_&, psdb& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::ssbn::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::V_th, V_th );
  def< double >( d, names::V_reset, V_reset );
  def< double >( d, names::t_ref, t_ref );
  def< double >( d, names::g_L, g_L );
  def< double >( d, names::E_L, E_L );
  def< double >( d, names::E_ex, E_ex );
  def< double >( d, names::E_in, E_in );
  def< double >( d, names::C_m, C_m );
  def< double >( d, names::tau_syn_ex, tau_synE );
  def< double >( d, names::tau_syn_in, tau_synI );
  def< double >( d, names::I_e, I_e );
  def< double >( d, names::spb, n_s);
}

void
nest::ssbn::Parameters_::set( const DictionaryDatum& d )
{
  // allow setting the membrane potential
  updateValue< double >( d, names::V_th, V_th );
  updateValue< double >( d, names::V_reset, V_reset );
  updateValue< double >( d, names::t_ref, t_ref );
  updateValue< double >( d, names::E_L, E_L );

  updateValue< double >( d, names::E_ex, E_ex );
  updateValue< double >( d, names::E_in, E_in );

  updateValue< double >( d, names::C_m, C_m );
  updateValue< double >( d, names::g_L, g_L );

  updateValue< double >( d, names::tau_syn_ex, tau_synE );
  updateValue< double >( d, names::tau_syn_in, tau_synI );

  updateValue< double >( d, names::I_e, I_e );
  updateValue< double >( d, names::spb, n_s );

  if ( V_reset >= V_th )
    throw BadProperty( "Reset potential must be smaller than threshold." );

  if ( C_m <= 0 )
    throw BadProperty( "Capacitance must be strictly positive." );

  if ( t_ref < 0 )
    throw BadProperty( "Refractory time cannot be negative." );

  if ( tau_synE <= 0 || tau_synI <= 0 )
    throw BadProperty( "All time constants must be strictly positive." );
}

void
nest::ssbn::State_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::V_m, y[ V_M ] ); // Membrane potential
}

void
nest::ssbn::State_::set( const DictionaryDatum& d,
  const Parameters_& )
{
  updateValue< double >( d, names::V_m, y[ V_M ] );
}


/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

nest::ssbn::psdb()
  : Archiving_Node()
  , P_()
  , S_( P_ )
  , B_( *this )
{
  recordablesMap_.create();
}

nest::ssbn::psdb( const psdb& n )
  : Archiving_Node( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

nest::ssbn::~psdb()
{
  // GSL structs may not have been allocated, so we need to protect destruction
  if ( B_.s_ )
    gsl_odeiv_step_free( B_.s_ );
  if ( B_.c_ )
    gsl_odeiv_control_free( B_.c_ );
  if ( B_.e_ )
    gsl_odeiv_evolve_free( B_.e_ );
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
nest::ssbn::init_state_( const Node& proto )
{
  const ssbn& pr = downcast< psdb >( proto );
  S_ = pr.S_;
  FILE *fd = fopen("/dev/urandom", "r");

  size_t seed;

  fread(&seed, sizeof(seed),1,fd);
  fclose(fd);
  //printf("%lu", seed);
  srand(seed);

}

void
nest::ssbn::init_buffers_()
{
  Archiving_Node::clear_history();

  B_.spike_exc_.clear(); // includes resize
  B_.spike_inh_.clear(); // includes resize
  B_.currents_.clear();  // includes resize

  B_.logger_.reset();

  B_.step_ = Time::get_resolution().get_ms();
  B_.IntegrationStep_ = B_.step_;

  if ( B_.s_ == 0 )
    B_.s_ =
      gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, State_::STATE_VEC_SIZE );
  else
    gsl_odeiv_step_reset( B_.s_ );

  if ( B_.c_ == 0 )
    B_.c_ = gsl_odeiv_control_y_new( 1e-3, 0.0 );
  else
    gsl_odeiv_control_init( B_.c_, 1e-3, 0.0, 1.0, 0.0 );

  if ( B_.e_ == 0 )
    B_.e_ = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  else
    gsl_odeiv_evolve_reset( B_.e_ );

  B_.sys_.function = ssbn_dynamics;
  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params = reinterpret_cast< void* >( this );

  B_.I_stim_ = 0.0;
}

void
nest::ssbn::calibrate()
{
  // ensures initialization in case mm connected after Simulate
  B_.logger_.init();

  V_.PSConInit_E = 1.0 * numerics::e / P_.tau_synE;
  V_.PSConInit_I = 1.0 * numerics::e / P_.tau_synI;
  V_.RefractoryCounts = Time( Time::ms( P_.t_ref ) ).get_steps();

  // since t_ref >= 0, this can only fail in error
  assert( V_.RefractoryCounts >= 0 );
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void
nest::ssbn::update( Time const& origin,
  const long from,
  const long to )
{
  assert(
    to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );
  //srand(time(NULL));
  for ( long lag = from; lag < to; ++lag )
  {

    double t = 0.0;

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals


	while ( t < B_.step_ )
	{
	  const int status = gsl_odeiv_evolve_apply( B_.e_,
	    B_.c_,
	    B_.s_,
	    &B_.sys_,             // system of ODE
	    &t,                   // from t
	    B_.step_,             // to t <= step
	    &B_.IntegrationStep_, // integration step size
	    S_.y );               // neuronal state


	  if ( status != GSL_SUCCESS )
	    throw GSLSolverFailure( get_name(), status );
	}

	if (S_.r)
	{
		--S_.r;
	      	S_.y[ State_::V_M ] = P_.V_reset; // clamp potential
	}

    // refractoriness and spike generation

     	/*		
    if ( S_.r )
    { // neuron is absolute refractory
      --S_.r;
      S_.y[ State_::V_M ] = P_.V_reset; // clamp potential
    }
	*/

	/*	
    if (S_.r ==0)
    {
	    S_.y[ State_::DG_EXC ] += B_.spike_exc_.get_value( lag ) * V_.PSConInit_E;
	    S_.y[ State_::DG_INH ] += B_.spike_inh_.get_value( lag ) * V_.PSConInit_I;

    }
    else	
      --S_.r;
	*/
    else if ( S_.y[ State_::V_M ] >= P_.V_th ) 
    //if ( S_.y[ State_::V_M ] >= P_.V_th ) 
    //if ( S_.y[ State_::V_M ] >= P_.V_th ) 
 // neuron is not absolute refractory
    {
      S_.y[ State_::V_M ] = P_.V_reset;
	
 	S_.rand_x = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
	/*
	if (S_.rand_x == 1.0 || S_.rand_x == 0.0)
	{
		S_.rand_x = static_cast<double>(std::rand())/static_cast<double>(RAND_MAX);
		//S_.rand_x = 0.4;
	}
	*/
 	//S_.rand_x = 1.0;
	//printf("rand_x:%4.2f",rand_x);
	if (S_.bm==0 && S_.rand_x < 1.0/P_.n_s)
	{
		S_.bm = P_.n_s -1;

      		S_.r = V_.RefractoryCounts;
	      // log spike with Archiving_Node
	      set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

	      SpikeEvent se;
	      kernel().event_delivery_manager.send( *this, se, lag );
        }
	/*
	else
	{
      		S_.r = V_.RefractoryCounts;
	}
	*/
    }

    else if (S_.bm>0 && S_.r == 0)
    {
        S_.r = V_.RefractoryCounts;                                                                 
        S_.y[ State_::V_M ] = P_.V_reset;                                                                         
        --S_.bm;                                                                                      
        set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );                                  
        SpikeEvent se;                                                                                
        kernel().event_delivery_manager.send( *this, se, lag );                                       
  
    }   

    
    // add incoming spikes
    S_.y[ State_::DG_EXC ] += B_.spike_exc_.get_value( lag ) * V_.PSConInit_E;
    S_.y[ State_::DG_INH ] += B_.spike_inh_.get_value( lag ) * V_.PSConInit_I;
    
    // set new input current
    B_.I_stim_ = B_.currents_.get_value( lag );

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
  }
}

void
nest::ssbn::handle( SpikeEvent& e )
{
  assert( e.get_delay() > 0 );

  if ( e.get_weight() > 0.0 )
    B_.spike_exc_.add_value( e.get_rel_delivery_steps(
                               kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  else
    B_.spike_inh_.add_value( e.get_rel_delivery_steps(
                               kernel().simulation_manager.get_slice_origin() ),
      -e.get_weight()
        * e.get_multiplicity() ); // ensure conductance is positive
}

void
nest::ssbn::handle( CurrentEvent& e )
{
  assert( e.get_delay() > 0 );

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(
    e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
    e.get_weight() * e.get_current() );
}

void
nest::ssbn::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

#endif // HAVE_GSL
