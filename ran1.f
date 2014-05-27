C-------------------------------------------------------------------
C     
C (PSEUDO) RANDOM NUMBER GENERATOR     
C ================================
C                            
C RETURNS A UNIFORM RANDOM DEVIATE BETWEEN 0.0 AND 1.0. 
C SET IDUM TO ANY NEGATIVE INTEGER VALUE TO INITIALIZE 
C OR REINITIALIZE THE SEQUENCE.
C                                                                      
C FROM "NUMERICAL RECIPES"                                            
C
C SAVE statement for IFF, R, IXn      B. Schmidt, May  31, 1992
C RAN1 returned as REAL*8             B. Schmidt, May  31, 1992
C
C-------------------------------------------------------------------
      REAL*8 FUNCTION RAN1(IDUM)    

      real*8 r(97),rm1,rm2                                                
      integer ix1,ix2,ix3
      integer ia1,ia2,ia3
      integer ic1,ic2,ic3
      integer m1,m2,m3
      integer idum,j,iff

      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247D-6)     
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773D-6)    
      PARAMETER (M3=243000,IA3=4561,IC3=51349)                    

      SAVE IFF, R, IX1,IX2,IX3    !?????????
      DATA IFF /0/                                               

      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN                           
         IFF=1                                      
         IX1=MOD(IC1-IDUM,M1)                      
         IX1=MOD(IA1*IX1+IC1,M1)                  
         IX2=MOD(IX1,M2)                         
         IX1=MOD(IA1*IX1+IC1,M1)                
         IX3=MOD(IX1,M3)                       
         DO 11 J=1,97                         
            IX1=MOD(IA1*IX1+IC1,M1)           
            IX2=MOD(IA2*IX2+IC2,M2)          
            R(J)=(dble(IX1)+dble(IX2)*RM2)*RM1
 11      CONTINUE                              
         IDUM=1                               
      ENDIF                                 
      IX1=MOD(IA1*IX1+IC1,M1)              
      IX2=MOD(IA2*IX2+IC2,M2)             
      IX3=MOD(IA3*IX3+IC3,M3)            
      J=1+(97*IX3)/M3                   
      IF(J.GT.97.OR.J.LT.1)STOP      
      RAN1=  R(J) 
      R(J)=(dble(IX1)+dble(IX2)*RM2)*RM1


      RETURN                             
      END                               

