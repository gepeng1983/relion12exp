/***************************************************************************
 *  *
 *   * Author: "Peng Ge"
 *    *  Univeristy of California Los Angeles
 *     *
 *      * This program is free software; you can redistribute it and/or modify
 *       * it under the terms of the GNU General Public License as published by
 *        * the Free Software Foundation; either version 2 of the License, or
 *         * (at your option) any later version.
 *          *
 *           * This program is distributed in the hope that it will be useful,
 *            * but WITHOUT ANY WARRANTY; without even the implied warranty of
 *             * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              * GNU General Public License for more details.
 *               *
 *                * This complete copyright notice must be included in any revised version of the
 *                 * source code. Additional authorship citations may be added, but existing
 *                  * author citations must be preserved.
 *                   * 
 *                    * This routine is written after IHRSR himpose by Edward Egelman of 
 *                     * Unversity of Virginia 
 *                      ***************************************************************************/

#ifndef HELIX_H_
#define HELIX_H_

#include "src/multidim_array.h"
#include "src/image.h"

#define dgr (2*PI/360.)

#define DEBUG

template <typename T>
class Heliciser 
{

public:

	T* dencyl;
	T* symcyl;

	T apix;

	T hrise, hturn, hinner, houter, htake;
	
	long int hsuper;
	long int irmin, irmax, nphi;
	
	T deltaz;

	long int nzlim;

	void himpose(MultidimArray<T> &img)
	{
		long int nrR=-irmin+irmax+1;
		dencyl=new T[nrR*nphi*ZSIZE(img)];

		memset(dencyl, 0, nrR*nphi*ZSIZE(img)*sizeof(T));
	
//		std::cout << "before work" << std::endl;
		work(img, dencyl);
//		std::cout << "after  work" << std::endl;

		

		symcyl=new T[nzlim*nphi*(irmax-irmin+1)+1024];

		memset(symcyl, 0, (nzlim*nphi*(irmax-irmin+1)+1024)*sizeof(T));

//		std::cout << "before sym " << std::endl;
		symmetry(img, dencyl, symcyl);
//		std::cout << "after sym  " << std::endl;

		delete dencyl;
		dencyl=NULL;

#ifdef DEBUG
		Image<double> im(nrR, nphi, nzlim);
		T* tmpdata=im.data.data;
		im.data.data=symcyl;
		im.write("symcyl.mrc");
		im.data.data=tmpdata;
#endif
	
		output(symcyl, img);
	
		delete symcyl;
		symcyl=NULL;
	
	}
//#define DEBUG

	void hsearch(MultidimArray<T> &img, T hriseinc, T hturninc, T * logx = NULL, T * logy = NULL)
	{
		long int nrR=-irmin+irmax+1;
		dencyl=new T[nrR*nphi*ZSIZE(img)];
		
		memset(dencyl, 0, nrR*nphi*ZSIZE(img)*sizeof(T));
		work(img, dencyl);

		
#ifdef DEBUG
		std::cout << "step 0 original hrise, hturn: " << hrise << ", " << hturn << std::endl;
#endif

		T x[21], y[21];
		
		//fix rise search turn
		for (int i=-10; i<=10; i++)
		{
			x[i+10]=hturn+i*hturninc;
			T residue = hresidue(img, dencyl, hrise, x[i+10]);
			y[i+10]=residue;
		}
		
		hturn = bottom(x, y);
		
		if (logx && logy)
			for (int i=0;i<21;i++)
			{
				logx[i]=x[i];
				logy[i]=y[i];
#ifdef DEBUG
				std::cout << x[i] << ", " << y[i] << std::endl;
#endif
			}
		
#ifdef DEBUG
		std::cout << "step 1 new hrise, hturn: " << hrise << ", " << hturn << std::endl;
#endif
		
		//fix turn search rise
		for (int i=-10; i<=10; i++)
		{
			x[i+10]=hrise+i*hriseinc;
			T residue = hresidue(img, dencyl, x[i+10], hturn);
			y[i+10]=residue;
		}
		
		hrise = bottom(x, y);
		
		if (logx && logy)
			for (int i=0;i<21;i++)
			{
				logx[i+21]=x[i];
				logy[i+21]=y[i];
#ifdef DEBUG
				std::cout << x[i] << ", " << y[i] << std::endl;
#endif
			}

#ifdef DEBUG
		std::cout << "step 2 new hrise, hturn: " << hrise << ", " << hturn << std::endl;
#endif
	
		//fix rise search turn
		for (int i=-10; i<=10; i++)
		{
			x[i+10]=hturn+i*hturninc;
			T residue = hresidue(img, dencyl, hrise, x[i+10]);
			y[i+10]=residue;
		}
		
		hturn = bottom(x, y);
		
		if (logx && logy)
			for (int i=0;i<21;i++)
			{
				logx[i+42]=x[i];
				logy[i+42]=y[i];
#ifdef DEBUG
				std::cout << x[i] << ", " << y[i] << std::endl;
#endif
			}

#ifdef DEBUG
		std::cout << "step 3 new hrise, hturn: " << hrise << ", " << hturn << std::endl;
#endif
		
		//fix turn search rise
		for (int i=-10; i<=10; i++)
		{
			x[i+10]=hrise+i*hriseinc;
			T residue = hresidue(img, dencyl, x[i+10], hturn);
			y[i+10]=residue;
		}
		
		hrise = bottom(x, y);
		
		if (logx && logy)
			for (int i=0;i<21;i++)
			{
				logx[i+63]=x[i];
				logy[i+63]=y[i];
#ifdef DEBUG
				std::cout << x[i] << ", " << y[i] << std::endl;
#endif
			}

#ifdef DEBUG
		std::cout << "step 4 new hrise, hturn: " << hrise << ", " << hturn << std::endl;
#endif
	
		delete dencyl;
		dencyl=NULL;
	}

	T bottom(T* x, T* y)
	{
		int besti=-10;
		T besty=y[0];
		for (int i=-10; i<=10; i++)
		{
			if (besty > y[i+10]) {
				besti=i;
				besty=y[i+10];
			}
		}
		return x[besti+10];
	}

	T hresidue(MultidimArray<T> &img, T* dencyl, T hrise1, T hturn1)
	{
		T phiinc=2.*PI/nphi;
		//variables not changing
		T subl=(htake* (((T)ZSIZE(img))*((T)hsuper))/((T)deltaz) );
		long int nsubl=(long int)floor(subl);
      		T ss=(T)hsuper;
		long int nrR=irmax-irmin+1;

		//variables changing
		T deltaz1=hrise1/apix*hsuper;
		T deltadeltaz1=deltaz1/((T)nzlim);
		T zbase1=(ZSIZE(img)*hsuper-nsubl*deltaz1)/2.;



		long int iz, iphi, ir, k;

//		std::cout << nphi << std::endl;

		T tsum=0.;

		for (iz=0;iz<nzlim;iz++)
			for (iphi=0;iphi<nphi;iphi++)
				for (ir=0; ir<irmax-irmin+1;ir++)
				{
//					std::cout << "ir, iphi, iz " << ir << "," << iphi << "," << iz << std::endl;

					T dsum=0.;
					T dsum2=0.;
	
					for (k=0;k<nsubl;k++)
					{
						T phip= iphi + (k*hturn1*dgr)/ phiinc;
//						std::cout << "k " << k << "phip " << phip << "moded " << fmod(phip,(T)nphi) << std::endl;
						phip=fmod(phip,(T)nphi);
						if (phip<0) phip+=(T)nphi;
						
						long int iphip=(long int)floor(phip);
						T dphip=phip-iphip;

						T z=(k*deltaz1 + iz*deltadeltaz1 + zbase1)/ss;
						long int izz=(long int)floor(z);
						T dizz=z-izz;
//						std::cout << "iphip " << iphip << " dphip " << dphip << " izz " << izz << " dizz " << dizz << std::endl;

/*						if (izz>ZSIZE(img)) {
							std::cerr << "Helix: iz =" << iz << std::endl;
							REPORT_ERROR("Helix: iz out of bound");
						}
*/
						T den1, den2, den3, den4;
						den1=dencyl[ir+iphip*nrR+izz*nrR*nphi];
						den2=dencyl[ir+iphip*nrR+(izz+1)*nrR*nphi];
						
						if (phip + 1.> (T) nphi) 
						{
							den3=dencyl[ir+izz*nrR*nphi];
							den4=dencyl[ir+(izz+1)*nrR*nphi];
						}
						else
						{
							den3=dencyl[ir+(iphip+1)*nrR+izz*nrR*nphi];
							den4=dencyl[ir+(iphip+1)*nrR+(izz+1)*nrR*nphi];
						}
						
						T d=(1.-dizz)*(1.-dphip)*den1 + dizz*(1.-dphip)*den2 + (1.-dizz)*dphip*den3 + dizz*dphip*den4;
						dsum+=d;
						dsum2+=d*d;
					}
					T sigma2=dsum2-(dsum*dsum/nsubl);
					
					tsum+=sigma2;
				}
		return tsum;

	}
//#undef DEBUG

	void work(MultidimArray<T> &img, T* dencyl)
	{
		T phiinc=2.*PI/nphi;

		long int k, j, i, iphi, ir;

		T* csphi=new T[nphi];
		T* snphi=new T[nphi];

/*	Use memset to set zeros
 *
 * 		for (k=0;k<ZSIZE(img);k++)
			for (j=0;j<nphi;j++)
				for (i=0;i<irmax-irmin+1;i++)
				{
					dencyl[i+j*nphi+k*nphi*ZSIZE(img)]=(T)0.;
				}
*/	
		for (i=0;i<nphi;i++)
		{
			T phi=i*phiinc;
			csphi[i]=cos(phi);
			snphi[i]=sin(phi);
		}

		T ss=(T) hsuper;

		long int nrR=irmax-irmin+1;

		img.setXmippOrigin();

//		std::cout << "nphi=" << nphi << std::endl;

		for (k=STARTINGZ(img);k<FINISHINGZ(img);k++)
			for (iphi=0;iphi<nphi;iphi++)
				for (ir=0;ir<irmax-irmin+1;ir++)
				{
					T r=(T)(ir+irmin);
					T x=r*csphi[iphi]/ss;
					T y=r*snphi[iphi]/ss;
	
					//std::cout << "Work: line 1 x, y, k=" << x << "," << y << "," << k << std::endl;
					//std::cout << "Work: line 1 ir, iphi=" << ir << "," << iphi << std::endl;
					//std::cout << "Work: index:" << ir+iphi*nrR+(k-STARTINGZ(img))*nrR*nphi << std::endl;
					dencyl[ir+iphi*nrR+(k-STARTINGZ(img))*nrR*nphi]=img.interpolatedElement3D(x,y,k);
				}
	
		delete csphi, snphi;
	}


	void symmetry(MultidimArray<T> &img, T* dencyl, T* symcyl)
	{
		T phiinc=2.*PI/nphi;
		T deltadeltaz=deltaz/((T)nzlim);
		T subl=(htake* (((T)ZSIZE(img))*((T)hsuper))/((T)deltaz) );
		long int nsubl=(long int)floor(subl);


		T zbase=(ZSIZE(img)*hsuper-nsubl*deltaz)/2.;

      		T ss=(T)hsuper;

		long int nrR=irmax-irmin+1;

		long int iz, iphi, ir, k;

//		std::cout << nphi << std::endl;
//
		std::cout << "deltaz " << deltaz << " deltadeltaz " << deltadeltaz << " zbase " << zbase << " ss " << ss << std::endl;

		for (iz=0;iz<nzlim;iz++)
			for (iphi=0;iphi<nphi;iphi++)
				for (ir=0; ir<irmax-irmin+1;ir++)
				{
					T dsum=0.;

//					std::cout << "ir, iphi, iz " << ir << "," << iphi << "," << iz << std::endl;
	
					for (k=0;k<nsubl;k++)
					{
						T phip= iphi + (k*hturn*dgr)/ phiinc;
//						std::cout << "k " << k << "phip " << phip << "moded " << fmod(phip,(T)nphi) << std::endl;
						phip=fmod(phip,(T)nphi);
						if (phip<0) phip+=(T)nphi;
						
						long int iphip=(long int)floor(phip);
						T dphip=phip-iphip;

						T z=(k*deltaz + iz*deltadeltaz + zbase)/ss;
						long int izz=(long int)floor(z);
						T dizz=z-izz;

//						if (dizz < 0 || dphip < 0) 
//						std::cout << "iphip " << iphip << " dphip " << dphip << " izz " << izz << " dizz " << dizz << std::endl;
//						std::cout << "iz " << iz << " k " << k << " z " << z << " izz " << izz << " dizz " << dizz << std::endl;

/*						if (izz>ZSIZE(img)) {
							std::cerr << "Helix: iz =" << iz << std::endl;
							REPORT_ERROR("Helix: iz out of bound");
						}
*/
						T den1, den2, den3, den4;
						den1=dencyl[ir+iphip*nrR+izz*nrR*nphi];
						den2=dencyl[ir+iphip*nrR+(izz+1)*nrR*nphi];
						
						if (phip + 1.> (T) nphi) 
						{
							den3=dencyl[ir+izz*nrR*nphi];
							den4=dencyl[ir+(izz+1)*nrR*nphi];
						}
						else
						{
							den3=dencyl[ir+(iphip+1)*nrR+izz*nrR*nphi];
							den4=dencyl[ir+(iphip+1)*nrR+(izz+1)*nrR*nphi];
						}
						
						dsum+=(1.-dizz)*(1.-dphip)*den1 + dizz*(1.-dphip)*den2 + (1.-dizz)*dphip*den3 + dizz*dphip*den4;
					}

					symcyl[ir+iphi*nrR+iz*nrR*nphi]=dsum/nsubl;
				}
	}
	
	void output(T* symcyl, MultidimArray<T> &out)
	{
		long int naniso=10*hsuper + nzlim*(2 + floor(YSIZE(out)*hsuper/deltaz));
		T* imageout=new T[XSIZE(out)*XSIZE(out)*naniso];

		memset(imageout, 0, XSIZE(out)*XSIZE(out)*naniso*sizeof(T));
		
		
		T phiinc=2.*PI/nphi;
      		T ss=(T)hsuper;
		T radlim1=irmin*irmin/ss/ss;
		T radlim2=irmax*irmax/ss/ss;
		T subl=(htake* (((T)ZSIZE(out))*((T)hsuper))/((T)deltaz) );
		long int nsubl = (long int) floor(subl);
		T zbase=((T)ZSIZE(out)*hsuper-nsubl*deltaz)/2.;

		

		long int xcen=XSIZE(out)/2;
		long int ycen=YSIZE(out)/2;

		long int nsubnew=(ZSIZE(out)*hsuper/deltaz) + 1;

		long int nsubbefore=(zbase)/deltaz+1;

		long int nrR=irmax-irmin+1;

		long int krep, kz, k, i, j;


		std::cout << "nsubnew " << nsubnew << " nsubbefore " << nsubbefore << " hsuper " << hsuper << " deltaz " << deltaz << " zbase " << zbase << std::endl;

		
//		std::cerr << "output line 1 " << xcen << " " << ycen << " " << hsuper << std::endl;
		
		for (krep=-nsubbefore; krep<=nsubnew-nsubbefore; krep++)
		{
//			std::cerr << "output line 2 krep=" << krep << std::endl;
			for (kz=0; kz<nzlim; kz++)
			{
				k=(krep+nsubbefore)*nzlim+kz;
//				std::cerr << "output line 3 krep, kz, k=" << krep << " " << kz << " " << k << std::endl;
				for (j=0; j<YSIZE(out); j++)
				{
					T y=j-ycen;
					T y2=y*y;
					for (i=0; i<XSIZE(out); i++)
					{

						T x=i-xcen;
						T r2=y2+x*x;
						
						if (r2> radlim2 || r2 < radlim1) continue;
			
						T phip = dgr*(1./dgr*atan2(y,x)-krep*hturn)/phiinc;
	
			//			std::cout << "phip, fmod(phip), r " << phip << "," << fmod(phip,nphi) << "," << sqrt(r2)*hsuper << std::endl;

						phip=fmod(phip,nphi);
						if (phip<0) phip+=(T)nphi;

			//			std::cout << "phip " << phip << std::endl;
						
						long int iphip=(long int)floor(phip);
						T dphip=phip-iphip;
						
						T r=sqrt(r2)*hsuper-irmin;
						long int ir=(long int)floor(r);
						T dir=r-ir;

//						std::cout << "iphip " << iphip << " dphip " << dphip << " ir " << ir << " dir " << dir << std::endl;
						

						T den1, den2, den3, den4;

						den1 = symcyl[ir+iphip*nrR+kz*nrR*nphi];
						den2 = symcyl[ir + 1+iphip*nrR+kz*nrR*nphi];

						if (phip + 1.> (T) nphi)
						{
							den3=symcyl[ir+kz*nrR*nphi];
							den4=symcyl[ir+1+kz*nrR*nphi];		
						}
						else
						{
							den3=symcyl[ir+(iphip+1)*nrR+kz*nrR*nphi];
							den4=symcyl[ir+1+(iphip+1)*nrR+kz*nrR*nphi];
						}
	
						T density=(1.-dir)*(1.-dphip)*den1 + dir*(1.-dphip)*den2 + (1.-dir)*dphip*den3 + dir*dphip*den4;
			//			std::cerr << "output line 4 i, j, k=" << i << "," << j << "," << k << std::endl;
//						if (kz==nzlim-1)
//							imageout[i+j*XSIZE(out)+k*XSIZE(out)*YSIZE(out)]=0;
//						else
						imageout[i+j*XSIZE(out)+k*XSIZE(out)*YSIZE(out)]=density;

					}
				}
			}
		}

#ifdef DEBUG
		Image<double> im(XSIZE(out),XSIZE(out),naniso);
		T* tmpdata=im.data.data;
		im.data.data=imageout;
		im.write("imageout.mrc");
		im.data.data=tmpdata;
#endif

		T zratio=((T)nzlim)/deltaz;
		T zbase2=nsubbefore*deltaz-zbase;

		out.initZeros();

//		std::cerr << "output line A" << std::endl;
		for (k=0;k<ZSIZE(out); k++)
		{
			T zwant=k*zratio*ss+zbase2*zratio;

			long int iz=(long int)floor(zwant);
			T zdif=zwant-iz;

			T omzdif=1.-zdif;

			long int izp1=iz+1;

			for (j=0; j< YSIZE(out); j++)
				for (i=0; i<XSIZE(out); i++)
				{
					T den1, den2, density;

					den1=imageout[i+j*XSIZE(out)+iz*XSIZE(out)*YSIZE(out)];	
					den2=imageout[i+j*XSIZE(out)+izp1*XSIZE(out)*YSIZE(out)];	
					
					density=omzdif*den1 + zdif*den2;
				
					DIRECT_NZYX_ELEM(out, 0, k, j, i) = density;
				}
		}
		delete imageout;
	}

	Heliciser(T newapix=1., T newhrise=1., T newhturn=1., T newhinner=-1., T newhouter=-1., long int newhsuper=1, T newhtake=0.5)
	{
		dencyl=NULL;
		symcyl=NULL;
		apix=newapix;
		hrise=newhrise;
		hturn=newhturn;
		hsuper=newhsuper;
		htake=newhtake;

		deltaz=hrise/apix*hsuper;
		nzlim=(long int)floor(deltaz + 2.0*hsuper);


		hinner=newhinner;
		houter=newhouter;
	
		if (newhinner > 0 || newhouter > 0) 
		{
			irmin=(long int)floor(hinner/apix*hsuper);
			irmax=(long int)floor(houter/apix*hsuper);
			nphi=(long int)floor(2.*PI*houter/apix*hsuper);
		}
		else
		{
			REPORT_ERROR("Helix: inner/outer radius not set");
		}
		
		
	}
	
	~Heliciser()
	{
		Clear();
	}
	
	void Clear()
	{
		if (dencyl) delete dencyl;
		if (symcyl) delete symcyl;
	}
};

#endif
