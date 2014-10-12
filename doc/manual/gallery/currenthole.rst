.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _gallery.mhdcurrenthole:


Current-Hole
^^^^^^^^^^^^

.. todo::

   needs to be validated

In this section we treat a special case of a 2D reduced MHD problem,
namely the Current-Hole problem. The general Current-Hole problem :cite:`Czarny2008` writes:

.. math::

  \left\{
  \begin{aligned}
  \partial_t \psi = (1+\epsilon x) [\psi,\phi] + \eta \left(  J - J_c\right),  \\
  \partial_t \omega = 2 \epsilon \frac{\partial \phi}{\partial y} \omega + (1+\epsilon x) [\omega, \phi] + \frac{1}{1+\epsilon x}[\psi,J]+ \nu \nabla^2 \omega, \\
  J = \Delta^{\star} \psi,  \\
  \Delta_{\perp} \phi = \omega .
  \end{aligned} 
  \right.

where :math:`[a,b] = \frac{\partial a}{\partial x_1} \frac{\partial b}{\partial x_2} - \frac{\partial a}{\partial x_2} \frac{\partial b}{\partial x_1}`, denotes the Poisson Bracket of the functions :math:`a,b`.

In the sequel, we consider only the planar cylindrical geometry (:math:`\epsilon=0`). The problem of the Current-Hole then writes: 

.. math::
  :label: current-hole-eq

  \left\{
  \begin{aligned}
  \partial_t \psi &= [\psi,\phi] + \eta \left(  J - J_c\right),
  \\
  \partial_t \omega &= [\omega, \phi] + [\psi,J]+ \nu \nabla^2 \omega,
  \\
  J &= \nabla^2 \psi,
  \\
  \nabla^2 \phi &= \omega.
  \end{aligned}  
  \right.

The current-hole problem, is subject to the initial conditions:

.. math::

  \left\{\begin{array}{c}
  J (0,\mathbf{x}) = J_c (\mathbf{x}) , ~~~ \mathbf{x} \in \Omega  \\
  \psi (0,\mathbf{x})  = \nabla^{-2} J_c , ~~~ \mathbf{x} \in \Omega  \\
  \phi (0,\mathbf{x})  = 0 , ~~~ \mathbf{x} \in \Omega  \\
  \omega (0,\mathbf{x}) = 0 , ~~~ \mathbf{x} \in \Omega  
  \end{array}\right.


and the boundary conditions:

.. math::

  \left\{\begin{array}{c}
  J (t,\mathbf{x}) = 0 , ~~~ \mathbf{x} \in \partial \Omega, ~t \in [0,T]  \\
  \psi (t,\mathbf{x})  = 0 , ~~~ \mathbf{x} \in \partial \Omega , ~t \in [0,T]  \\
  \phi (t,\mathbf{x})  = 0 , ~~~ \mathbf{x} \in  \partial \Omega, ~t \in [0,T]  \\
  \omega (t,\mathbf{x}) = 0 , ~~~ \mathbf{x} \in \partial \Omega, ~t \in [0,T]    
  \end{array}\right.


Before going to numerical simulations, let us recall some properties for the current-hole problem. We refer to :cite:`Deriaz2011`,:cite:`depres_sart_mhd` for proofs. 

\begin{proposition}[Conservation]
If :math:`\nu=\eta=0` then all regular solutions of (Eq. \ref{current_hole_eq}), verify:

* *Energy conservation*:

.. math::

  \frac{d}{dt}\int_{\Omega} | \nabla \psi|^2 + | \nabla \phi|^2 ~d\Omega &= 0,
 
* *Magnetic Helicity conservation* :

.. math::

  \frac{d}{dt}\int_{\Omega} \psi ~d\Omega &= 0,

* *Cross-Helicity conservation* :

.. math::  
  
  \frac{d}{dt}\int_{\Omega} \nabla \psi \cdot \nabla \phi ~d\Omega = - \frac{d}{dt}\int_{\Omega} \psi \phi ~d\Omega = 0.

\end{proposition}

As mentioned in :cite:`Deriaz2011`, when :math:`J_c=0` and :math:`\nu>0` and :math:`\eta>0`, we have:

.. math::

  \frac{d}{dt}\int_{\Omega} | \nabla \psi|^2 + | \nabla \phi|^2 ~d\Omega & \leq 0.

In the sequel, we present the time scheme and variational formulation that we used.


Time scheme
___________


We use a semi-implicit time scheme:

.. math::
  :label: mhd-currenthole-eq-w-LF

  \frac{\omega^{n+1} - \omega^{n}}{\Delta t} &=  [ \omega^{n} , \phi^{n} ] + [ \psi^{n} , J^{n} ] + \nu \nabla^2 \omega^{n+1}, 

.. math::  
  :label: mhd-currenthole-eq-phi-LF

  \nabla^2 \phi^{n+1} &=\omega^{n+1}, 

.. math::
  :label: mhd-currenthole-eq-psi-LF

  \frac{\psi^{n+1} - \psi^{n}}{\Delta t} &= [\psi^{n},\phi^{n+1}] + \eta  \nabla^2 \psi^{n+1} -  \eta  J_c, 

.. math::  
  :label: mhd-currenthole-eq-J-LF 

  J^{n+1} &=   \nabla^2 \psi^{n+1}. 

This time scheme has the advantage of being very simple, however it is of order 1.


Variational formulation
_______________________


Now let us introduce the discrete space, 

.. math::

  \mathcal{V}_h^0 &= \mathbf{span} \{ \varphi_b,~~~b \in \Lambda^0 \}.

where, the functions :math:`\varphi_b` can be *B-splines* or more generally *NURBS*.

Multiplying (Eq. \ref{mhd_currenthole_eq_w_LF}) by :math:`\varphi_b` and taking the integral over the whole domain, we get:

.. math::

  \frac{\int_{\Omega} \omega^{n+1} \varphi_b - \int_{\Omega} \omega^{n} \varphi_b}{\Delta t} &= 
  \int_{\Omega} [ \omega^{n} , \phi^{n} ]\varphi_b 
  + \int_{\Omega} [ \psi^{n} , J^{n} ]\varphi_b 
  + \nu \int_{\Omega} \nabla^2 \omega^{n+1} \varphi_b.

Now using Green's Formula and the boundary condition, we get:

.. math::

  \frac{\int_{\Omega} \omega^{n+1} \varphi_b - \int_{\Omega} \omega^{n} \varphi_b}{\Delta t} &= 
  \int_{\Omega} [ \omega^{n} , \phi^{n} ]\varphi_b
  + \int_{\Omega} [ \psi^{n} , J^{n} ]\varphi_b - \nu \int_{\Omega} \nabla \omega^{n+1} \cdot \nabla \varphi_b,

thus,

.. math::

  \int_{\Omega} \omega^{n+1} \varphi_b + \nu \Delta t \int_{\Omega} \nabla \omega^{n+1} \cdot \nabla \varphi_b  &= 
  \int_{\Omega} \omega^{n} \varphi_b 
  + \Delta t \int_{\Omega} [ \omega^{n} , \phi^{n} ]\varphi_b 
  +  \Delta t \int_{\Omega} [ \psi^{n} , J^{n} ]\varphi_b.

Let us consider :math:`\omega_h \in \mathcal{V}_h^{0}` an approximation of :math:`\omega`. We can expand :math:`\omega_h` over the basis of :math:`\mathcal{V}_h^0`:

.. math::

  \omega_h &= \sum_{b^{\prime} \in \Lambda^0} [\omega]^{b^{\prime}} \varphi_{b^{\prime}}.

Therefore, the previous equation leads to the linear system:

.. math::

  A_{\nu}^0 [\omega^{n+1}] = M^0[\omega^{n}] + \Delta t \mathcal{C}_{\mathcal{V}_h^0}[J^{n},\psi^{n}] + \Delta t \mathcal{C}_{\mathcal{V}_h^0}[\phi^{n}, \omega^{n}]

where we have introduced the matrices,

.. math::

  M^0 &= ( \int_{\Omega} \varphi_b \varphi_{b^{\prime}})_{b,b^{\prime} \in \Lambda^0},
  \\
  S^0 &= ( \int_{\Omega} \nabla \varphi_b \cdot \nabla  \varphi_{b^{\prime}})_{b,b^{\prime} \in \Lambda^0},
  \\
  A_{\nu}^0 &= M^0 + \Delta t \nu S^0

and :math:`\mathcal{C}_{\mathcal{V}_h^0}[u,v]` is the :math:`L^2` contribution over :math:`\mathcal{V}_h^0` of the Poisson's Bracket :math:`[u,v]`, \textit{i.e.} a column vector where the element of each line :math:`b` is :math:`\int_{\Omega} [ u,v ]\varphi_b`.

For (Eq. \ref{mhd_currenthole_eq_phi_LF}), we have for any :math:`\varphi_b \in \mathcal{V}_h^0`:

.. math::

  \int_{\Omega} \nabla \phi^{n+1} \cdot \nabla \varphi_b &= - \int_{\Omega} \omega^{n+1} \varphi_b

which leads to the linear system:

.. math::

  S^0 [\phi^{n+1}] = - M^0 [\omega^{n+1}]~.

Now, let's go back to (Eq. \ref{mhd_currenthole_eq_psi_LF}), we have for any :math:`\varphi_b \in \mathcal{V}_h^0`:

.. math::

  \frac{\int_{\Omega} \psi^{n+1} \varphi_b - \int_{\Omega} \psi^{n}\varphi_b}{\Delta t} &= 
  \int_{\Omega} [\psi^{n},\phi^{n+1}]\varphi_b 
  + \eta \int_{\Omega} \nabla^2 \psi^{n+1}\varphi_b 
  -  \eta \int_{\Omega} J_c \varphi_b

using Green's Formula and the boundary conditions we get:

.. math::

  \frac{\int_{\Omega} \psi^{n+1} \varphi_b - \int_{\Omega} \psi^{n}\varphi_b}{\Delta t} &= 
  \int_{\Omega} [\psi^{n},\phi^{n+1}]\varphi_b 
  - \eta \int_{\Omega} \nabla \psi^{n+1} \cdot \nabla  \varphi_b 
  - \eta \int_{\Omega} J_c \varphi_b~.

Let us consider :math:`\psi_h \in \mathcal{V}_h^{0}` an approximation of :math:`\psi`. We can expand :math:`\psi_h` over the basis of :math:`\mathcal{V}_h^0`:

.. math::

  \psi_h &= \sum_{b^{\prime} \in \Lambda^0} [\psi]^{b^{\prime}} \varphi_{b^{\prime}}

which leads to:

.. math::

  A_{\eta}^0 [\psi^{n+1}] &= 
  M^0[\psi^{n}] 
  + \Delta t \mathcal{C}_{\mathcal{V}_h^0}[ \psi^{n} , \phi^{n+1} ] 
  - \eta \Delta t M^0 [J_c]

with,

.. math::

  A_{\eta}^0 = M^0 + \Delta t \eta S^0~.

Finally, (Eq. \ref{mhd_currenthole_eq_J_LF}) leads to:

.. math::

  M^0 [J^{n+1}] = - S^0 [\psi^{n+1}]~.

Finally, the discretization of the incompressible \textit{MHD} writes:

.. math::
  :label: current-hole-eq1

  A_{\nu}^0 [\omega^{n+1}] &= M^0[\omega^{n}] + \Delta t \mathcal{C}_{\mathcal{V}_h^0}[J^{n},\psi^{n}] + \Delta t \mathcal{C}_{\mathcal{V}_h^0}[\phi^{n}, \omega^{n}],

.. math::
  :label: current-hole-eq2

  S^0 [\phi^{n+1}] &= - M^0 [\omega^{n+1}],

.. math::
  :label: current-hole-eq3

  A_{\eta}^0 [\psi^{n+1}] &= 
  M^0[\psi^{n}] 
  + \Delta t \mathcal{C}_{\mathcal{V}_h^0}[ \psi^{n} , \phi^{n+1} ] 
  - \eta \Delta t M^0 [J_c],

.. math::
  :label: current-hole-eq4

  M^0 [J^{n+1}] &= - S^0 [\psi^{n+1}].


As for the previous examples, we'll need to define for each equation two operators, one for the explicit term (for **un**, *i.e.* :math:`u_{n}`) and the other for
the implicit one (for **unew**, *i.e.* :math:`u_{n+1}`)

Hence for the equation :eq:`current-hole-eq1` we'll need::

  Ew = basicPDE(geometry=geo, testcase=tc_Ew)
  Iw = basicPDE(geometry=geo, testcase=tc_Iw)

  wn   = Ew.rhs
  wnew = Iw.unknwon

The same thing for :eq:`current-hole-eq2`, we get::   

  Ephi = basicPDE(geometry=geo, testcase=tc_Ephi)
  Iphi = basicPDE(geometry=geo, testcase=tc_Iphi)

  phin   = Ephi.rhs
  phinew = Iphi.unknwon

The same thing for :eq:`current-hole-eq3`, we get::   

  Epsi = basicPDE(geometry=geo, testcase=tc_Epsi)
  Ipsi = basicPDE(geometry=geo, testcase=tc_Ipsi)

  psin   = Epsi.rhs
  psinew = Ipsi.unknwon

And finally, for :eq:`current-hole-eq4`, we get::   

  Ej = basicPDE(geometry=geo, testcase=tc_Ej)
  Ij = basicPDE(geometry=geo, testcase=tc_Ij)

  jn   = Ej.rhs
  jnew = Ij.unknwon

Now let us define the dictionaries::

  # Ew dictionary
  tc_Ew = {}
  tc_Ew['AllDirichlet'] = True
  tc_Ew['b'] = lambda x,y : [1.]

  # Iw dictionary
  tc_Iw = {}
  tc_Iw['AllDirichlet'] = True
  tc_Iw['A'] = lambda x,y : [dt*nu, 0., 0., dt*nu]
  tc_Iw['b'] = lambda x,y : [1.]

  # Ephi dictionary
  tc_Ephi = {}
  tc_Ephi['AllDirichlet'] = True
  tc_Ephi['b'] = lambda x,y : [-1.]

  # Iphi dictionary
  tc_Iphi = {}
  tc_Iphi['AllDirichlet'] = True
  tc_Iphi['A'] = lambda x,y : [1., 0., 0., 1.]

  # Epsi dictionary
  tc_Epsi = {}
  tc_Epsi['AllDirichlet'] = True
  tc_Epsi['b'] = lambda x,y : [1.]

  # Ipsi dictionary
  tc_Ipsi = {}
  tc_Ipsi['AllDirichlet'] = True
  tc_Ipsi['A'] = lambda x,y : [dt*eta, 0., 0., dt*eta]
  tc_Ipsi['b'] = lambda x,y : [1.]

  # Ej dictionary
  tc_Ej = {}
  tc_Ej['AllDirichlet'] = True
  tc_Ej['A'] = lambda x,y : [-1., 0., 0., -1.]

  # Ij dictionary
  tc_Ij = {}
  tc_Ij['AllDirichlet'] = True
  tc_Ij['b'] = lambda x,y : [1.]

Now, we need to redefine the right hand side functions, in order to take into account the non linearity::
  
  # redefine the right hand side function for Ew
  def Fw(x,y):
        dxw  , dyw   = wn.grad(patch_id=0)
        dxj  , dyj   = jn.grad(patch_id=0)
        dxpsi, dypsi = psin.grad(patch_id=0)
        dxphi, dyphi = phin.grad(patch_id=0)

        v  = ( dxj * dypsi - dyj * dxpsi )
        v += ( dxphi * dyw - dyphi * dxw )
        v *= dt

        return [v]

  Ew.rhs.set_func(Fw)        
  
  # redefine the right hand side function for Epsi
  def Fpsi(x,y):
        dxpsi, dypsi = psin.grad(patch_id=0)
        dxphi, dyphi = phinew.grad(patch_id=0)

        v  = ( dxpsi * dyphi - dypsi * dxphi )
        v *= dt
        v -= eta * dt * jc.get() 

        return [v]

  Epsi.rhs.set_func(Fpsi) 

Finally, a single iteration of the Current-Hole problem can take the following form::

  # update w
  Ew.update()
  rhs = Ew.dot(wn) + Ew.rhs
  Iw.solve(rhs)

  # update phi
  rhs = Ephi.dot(wnew)
  Iphi.solve(rhs)

  # update psi
  Epsi.update()
  rhs = Epsi.dot(psin) + Epsi.rhs
  Ipsi.solve(rhs)

  # update j
  rhs = Ej.dot(psinew)
  Ij.solve(rhs)

.. Local Variables:
.. mode: rst
.. End:
