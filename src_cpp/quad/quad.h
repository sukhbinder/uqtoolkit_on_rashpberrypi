/* =====================================================================================
                      The UQ Toolkit (UQTk) version 2.1.1
                     Copyright (2013) Sandia Corporation
                      http://www.sandia.gov/UQToolkit/

     Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
     with Sandia Corporation, the U.S. Government retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is free software: you can redistribute it and/or modify
     it under the terms of the GNU Lesser General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
 
     You should have received a copy of the GNU Lesser General Public License
     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
 
     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
#ifndef QUAD_H_SEEN
#define QUAD_H_SEEN

#include "Array1D.h"
#include "Array2D.h"

/// \class	Quad
/// \brief	Generates quadrature rules
/// \note       Besides quadrature rules corresponding to PC bases, Clenshaw-Curtis(CC) and Newton-Cotes(NC) 
/// as well as ther Open (with no endpoints) versions are implemented.

class Quad {
public:
  
  Quad(Array1D<string>& grid_types, char *fs_type, Array1D<int>& param,Array1D<double>& alphas, Array1D<double>& bettas);
  /// \brief Constructor: initializes the rule type, sparseness type, dimensionality, level or grid size parameter, and two optional parameters for quadrature rule
  Quad(char *grid_type, char *fs_type,int ndim,int param,double alpha=0.0, double betta=1.0);
  /// \brief Constructor: empty
  Quad() {};
  /// \brief Destructor
  ~Quad() {}
  void init();

  /// \brief Set the parameter alpha
  void SetAlpha(double alpha){this->alpha_=alpha; }
  /// \brief Set the parameter beta
  void SetBeta(double betta){this->beta_=betta; }

  /// \brief Set the domain endpoints (for compact support domains)
  void SetDomain(Array1D<double>& aa, Array1D<double>& bb);
  /// \brief Set the domain endpoint (for semi-infinite domains)
  void SetDomain(Array1D<double>& aa);
  /// \brief Get the domain endpoints (for compact support domains)
  void GetDomain(Array1D<double>& aa, Array1D<double>& bb) const {aa=aa_; bb=bb_;}
  /// \brief Get the domain endpoint (for semi-infinite domains)
  void GetDomain(Array1D<double>& aa) const {aa=aa_;}
  
  /// \brief Set the rule externally (only quadrature points and weights)
  void SetRule(Array2D<double>& q, Array1D<double>& w);
  /// \brief Set the rule externally (quadrature points, weights and indices)
  /// Dummy function for backward compatibility
  void SetRule(Array2D<double>& q, Array1D<double>& w, Array2D<int>& ind){this->SetRule(q,w);}
  /// \brief Set the rule externally (quadrature points, weights, indices, and the level)
  //void SetRule(Array2D<double>& q, Array1D<double>& w, Array2D<int>& ind, int param);
  
  /// \brief Set the rule (the function that builds quadrature points/weights/indices)
  void SetRule();

  /// \brief Get the quadrature rule
  void GetRule(Array2D<double>& q, Array1D<double>& w);
  /// \brief Get the quadrature rule with indexing
    /// Dummy function for backward compatibility
    void GetRule(Array2D<double>& q, Array1D<double>& w, Array2D<int>& ind){this->GetRule(q,w);}
    
  /// \brief Compute the current points' indices at the next level
  //void ComputeCurrIndicesAtNextLevel(Array2D<int>& indAtNextLevel);

  /// \brief Externally set quadrature points
  void SetQdpts(Array2D<double>& q){  rule_.qdpts=q; return;}
  /// \brief Externally set the weights
  void SetWghts(Array1D<double>& w){  rule_.wghts=w; return;}
  /// \brief Externally set the indices
  // void SetIndices(Array2D<int>& ind){  rule_.indices=ind; return;}
  
  /// \brief Get quadrature points
  void GetQdpts(Array2D<double>& q){ q=rule_.qdpts;  return;}
  /// \brief Get the weights
  void GetWghts(Array1D<double>& w){  w=rule_.wghts; return;}
  /// \brief Get the indices
  //void GetIndices(Array2D<int>& ind){  ind=rule_.indices; return;}

  /// \brief Set the level parameter
  void SetLevel(int param) {nlevel_=param; return;}
    
  /// \brief Compute the indices of the next-level points
  void nextLevel();
  
  /// \brief Compress the rule, i.e. merge repeating points
  void compressRule();

  /// \brief Get the number of quadrature points
  int GetNQ() {return rule_.qdpts.XSize(); }
  
  /// \brief Get the indexing flag (whether indexing is turned on or off)
  // bool GetIndexing() {return indexing_;}
  
 private:
  
  /// \brief Dummy copy constructor, which should not be used as it is currently not well defined
  Quad(const Quad &) {};
          
  /// \brief The left endpoints of the domain
  Array1D<double> aa_;
  /// \brief the right endpoints of the domain
  Array1D<double> bb_;
  
  /// \brief The first parameter of the rule, if any
  double alpha_;
  /// \brief The second parameter of the rule, if any
  double beta_;
 /// \brief The first parameter of the rule, if any
  Array1D<double> alphas_;
  /// \brief The second parameter of the rule, if any
  Array1D<double> betas_;

  /// \brief Rule structure that stores quadrature points, weights and indices
  typedef struct 
  {
    /// \brief Quadrature points
    Array2D<double> qdpts;
    /// \brief Quadrature weights
    Array1D<double> wghts;
    /// \brief Quadrature indices (useful only for nested rules)
    //Array2D<int> indices;
  } QuadRule;
  
  /// \brief The quadrature rule structure
  QuadRule rule_;

  /// \brief The dimensionality
  int ndim_;
  
  /// \brief The current level, working variable for hierarchical construction
  int nlevel_;
  
  /// \brief The level for sparse rules, or the number of grid points per dim 
  /// for full product rules
  int maxlevel_;
  Array1D<int> param_;
  
  /// \brief Multiply two rules (full tensor product)
  void MultiplyTwoRules(QuadRule *rule1,QuadRule *rule2,QuadRule *rule_prod);
  /// \brief Multiply many rules (full tensor product)
  void MultiplyManyRules(int nrules, QuadRule *rules, QuadRule *rule_prod);
  /// \brief Add two rules
  void AddTwoRules(QuadRule *rule1,QuadRule *rule2,QuadRule *rule_sum);
  /// \brief Subtract two rules
  void SubtractTwoRules(QuadRule *rule1,QuadRule *rule2,QuadRule *rule_sum);
  
  /// \brief Compute 1D rules
  void create1DRule(string gridtype,Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b);

  /// \brief Clenshaw-Curtis (includes the endpoints)
  /// \note Heavily adopted from http://people.sc.fsu.edu/~jburkardt/cpp_src/sparse_grid_cc/sparse_grid_cc.html (distributed under LGPL)
  void create1DRule_CC(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b);
  /// \brief Legendre-Uniform
  void create1DRule_LU(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b);
  /// \brief Gauss-Hermite
  void create1DRule_HG(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr);
  /// \brief Newton-Cotes (i.e. equispaced, includes the endpoints)
  void create1DRule_NC(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b);
  /// \brief Newton-Cotes open (i.e. excludes the endpoints)
  void create1DRule_NCO(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b);
  /// \brief Clenshaw-Curtis open (i.e. excludes the endpoints)
  /// \note Heavily adopted from http://people.sc.fsu.edu/~jburkardt/cpp_src/sparse_grid_cc/sparse_grid_cc.html (distributed under LGPL)
  void create1DRule_CCO(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b);
  /// \brief Jacobi-Beta
  void create1DRule_JB(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b);
  /// \brief Gamma-Laguerre
  void create1DRule_GLG(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr);
  /// \brief Stieltjes-Wigert 
  void create1DRule_SW(Array1D<double>& qdpts,Array1D<double>& wghts,int ngr);
  /// \brief Custom rule given the recursive coefficients of the corresponding orthogonal polynomials
  /// \todo Recursive coefficients are given in a file 'ab.dat'; will need to make this more friendly
  void create1DRule_pdf(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b);
  /// \brief Gauss-Patterson starting with Legendre-Uniform 3
  void create1DRule_GP3(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b);
  
  /// \brief Auxilliary function: get the level of the multi-index
  void getMultiIndexLevel(Array2D<int>& multiIndexLevel, int level, int ndim);

  /// \brief Compute the number of grid points given the level and the growth rule
  //int grid(int il,int growthrule);
  //int grid(int il);

  /// \brief Indexing flag indicates whether indexing is performed or not 
  /// it is useful for nested sparse rules     
  //bool indexing_;
	
  /// \brief Growth rule: exponential(0) or linear(1)
  int growth_rule_;
  /// \brief Growth rules: exponential(0) or linear(1)
  Array1D<int> growth_rules_;
  
  /// \brief Grid type: 'CC','CCO','NC','NCO','LU', 'HG', 'JB', 'GLG', 'SW', 'pdf', or GP3
  string grid_type_;
  /// \brief Vector of grid types: 'CC','CCO','NC','NCO','LU', 'HG', 'JB', 'GLG', 'SW', 'pdf', or 'GP3'
  Array1D<string> grid_types_;
  
  /// \brief Sparseness type (full or sparse)
  string fs_type_;
  
}; 
#endif /* QUAD_H_SEEN */
