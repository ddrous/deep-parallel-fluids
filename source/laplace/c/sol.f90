      subroutine sol  (vkgs,vkgd,vkgi,vfg,kld,vu,neq,mp,ifac,isol, &
       nsym,energ,ier)

!   resolution d'un systeme lineaire symetrique ou non. la matrice est
!   stockee par ligne de ciel,en memoire dans les tables vkgs,vkgd,vkgi
!   il faut avoir reserve assez de memoire pour vkgs et vkgi
!   vkgi n'est pas utilise si la matrice est symetrique
!
!       entrees
!          vkgs,vkgd,vkgi    matrice du systeme : parties superieure,
!                            diagonale, inferieure (double precision)
!          vfg               second membre
!          kld               pointeurs vers les hauts de colonne
!          vu                vecteur solution (qui peut etre vfg)
!          neq               nombre d'equations
!          mp                unite logique d'impression (6 pour affichage sur ecran)
!          ifac              si ifac.eq.1 triangularisation de
!                            la matrice
!          isol              si isol.eq.1 calcul de la solution a
!                            partir de la matrice triangularisee
!          nsym              indice de probleme non symetrique (0:sym   1:non sym)
!       sorties
!          vkgs,vkgd,vkgi    matrice triangularisee (si ifac.eq.1)
!          vfg               solution (si isol.eq.1)
!          energ             energie du systeme (si nsym.eq.0)
!          ier               mis a 1 si pivot nul rencontre
!
!=========================== debut des declarations ====================
      implicit real*8 (a-h,o,q-z)
      implicit integer (p)
      dimension vkgs(*),vkgd(neq),vkgi(*),vfg(neq),kld(neq+1),vu(neq)
      data vzero/0.d0/ 
!=========================== debut du code executable ==================
!
!-------  traitement
!
      ik=1
      if(vkgd(1).eq.vzero) goto 800
      energ=vzero
      ier=0
      if (isol.eq.1) then
        do i = 1, neq
          vu(i) = vfg(i)
        end do
      end if
!
!-------  pour chaque colonne ik a modifier
!
      jhk=1
      do 100 ik=2,neq
!
!-------  pointeur du haut de la colonne suivante ik+1
!
      jhk1=kld(ik+1)
!
!-------  hauteur de la colonne ik (hors termes superieur et diagonal)
!
      lhk=jhk1-jhk
      lhk1=lhk-1
!
!-------  ligne du premier terme a modifier dans la colonne ik
!
      imin=ik-lhk1
      imin1=imin-1
!
!-------  ligne du dernier terme a modifier dans la colonne ik
!
      imax=ik-1
      if(lhk1.lt.0) goto 100
      if(ifac.ne.1) goto 90
      if(nsym.eq.1) vkgi(jhk)=vkgi(jhk)/vkgd(imin1)
      if(lhk1.eq.0) goto 40
!
!-------  modifier les termes non diagonaux de la colonne ik
!
      jck=jhk+1
      jhj=kld(imin)
!
!-------  pour chaque terme place en jck, correspondant a la colonne ij
!
      do 30 ij=imin,imax
      jhj1=kld(ij+1)
!
!-------  nombre de termes modificatifs du terme place en jck
!
      ic=min0(jck-jhk,jhj1-jhj)
      if(ic.le.0.and.nsym.eq.0) goto 20
      c1=vzero
      if(ic.le.0) goto 17
      j1=jhj1-ic
      j2=jck-ic
      if(nsym.eq.1) goto 15
      vkgs(jck)=vkgs(jck)-scal(vkgs(j1),vkgs(j2),ic)
      goto 20
15    vkgs(jck)=vkgs(jck)-scal(vkgi(j1),vkgs(j2),ic)
      c1=scal(vkgs(j1),vkgi(j2),ic)
17    vkgi(jck)=(vkgi(jck)-c1)/vkgd(ij)
20    jck=jck+1
30    jhj=jhj1
!
!-------  modifier le terme diagonal
!
40    jck=jhk
      cdiag=vzero
      do 70 ij=imin1,imax
      c1=vkgs(jck)
      if(nsym.eq.1) goto 50
      c2=c1/vkgd(ij)
      vkgs(jck)=c2
      goto 60
50    c2=vkgi(jck)
60    cdiag=cdiag+c1*c2
70    jck=jck+1
      vkgd(ik)=vkgd(ik)-cdiag
      if (vkgd(ik).eq.0.) goto 800
!
!-------  resolution du systeme triangulaire inferieur
!
90    if(isol.ne.1) goto 100
      if(nsym.ne.1) vu(ik)=vfg(ik)-scal(vkgs(jhk),vu(imin1),lhk)
      if(nsym.eq.1) vu(ik)=vfg(ik)-scal(vkgi(jhk),vu(imin1),lhk)
100   jhk=jhk1
      if(isol.ne.1) goto 9999
!
!-------  resolution du systeme diagonal :
!
      if(nsym.eq.1) goto 120
      do 110 ik=1,neq
      c1=vkgd(ik)
      if (c1.eq.vzero) goto 800
      c2=vu(ik)/c1
      vu(ik)=c2
110   energ=energ+c1*c2*c2
!
!-------  resolution du systeme triangulaire superieur
!
120   ik=neq+1
      jhk1=kld(ik)
130   ik=ik-1
      if(nsym.eq.1) vu(ik)=vu(ik)/vkgd(ik)
      if(ik.eq.1) goto 9999
      c1=vu(ik)
      jhk=kld(ik)
      jbk=jhk1-1
      if(jhk.gt.jbk)goto 150
      ij=ik-jbk+jhk-1
      do 140 jck=jhk,jbk
      vu(ij)=vu(ij)-vkgs(jck)*c1
140   ij=ij+1
150   jhk1=jhk
      goto 130
!
!-------  erreurs
!
800   if (mp.ne.0) write(mp,8000) ik
8000  format(' * sol pivot nul equation',i5)
      ier=1
      goto 9999
!
!-------  fin
!
9999  continue
      return
!===========================   fin du module sol    ==================
      end

      function scal(x,y,n)
!=======================================================================
! calcul du produit scalaire
!=======================================================================
      implicit real*8 (a-h,o-z)
      dimension x(1),y(1)
      data zero/0.d0/
!-----------------------------------------------------------------------
      scal=zero
      do i=1,n
      scal=scal+x(i)*y(i)
      enddo
      return
      end
