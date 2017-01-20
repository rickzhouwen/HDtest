        subroutine code1(n,a,c)
        implicit double precision (a-h,p-z)
        integer n
        
        double precision a(n,n)
        
        c=0.0
        do 100 i=1,n-1 
         do 200 j=i+1,n
     
            c=c+a(i,j)**2.d0   
          
200    continue        
100    continue        
       return 
       end
        

        subroutine code2(n,A,C)
        implicit double precision (a-h,p-z)
        integer n
        
        double precision A(n,n)
        
        C=0.0
        do 100 i=1,n 
         do 200 j=1,(n-1)
          if (j .NE. i) then        
            do 300 k=j+1,n   
               if (k .NE. i) then
                 C=C+A(i,j)*A(i,k)
               end if
                
300    continue
            end if
200    continue        
100    continue        
        
       end
       


        subroutine code3(n,A,C)
        implicit double precision (a-h,p-z)
        integer n
        
        double precision A(n,n)
        
        C=0.0
        do 100 i=1,n-3 
         do 200 j=i+1,n
    
            do 300 k=i+1,n-1
              if (k .NE. j) then
                do 400 l=k+1,n                   
           if ( l .NE. j) then
                 C=C+A(i,j)*A(k,l)                 
                 end if
                                              
400    continue 
            end if               
300    continue
200    continue        
100    continue        
        
       end




        subroutine code4(n,m,A,C)
        implicit double precision (a-h,p-z)
        integer n,m
        
        double precision A(n,m)
        
        C=0.0
        do 100 i=1,n 
         do 200 j=1,m
   
            C=C+A(i,j)**2.d0 

200    continue        
100    continue        
        
       end



        subroutine code5(n,m,A,C)
        implicit double precision (a-h,p-z)
        integer n,m
        
        double precision A(n,m)
        
        C=0.0
        do 100 i=1,n 
         do 200 j=1,m-1     
            do 300 k=j+1,m   
                 C=C+A(i,j)*A(i,k)               
300    continue
200    continue        
100    continue        
        
       end



        subroutine code6(n,m,A,C)
        implicit double precision (a-h,p-z)
        integer n,m
        
        double precision A(n,m)
        
        C=0.0
        do 100 i=1,n-1 
         do 200 j=1,m
    
            do 300 k=i+1,n
                do 400 l=1,m   
               if ( l .NE. j) then
                 C=C+A(i,j)*A(k,l)
               end if                                
400    continue              
300    continue
200    continue        
100    continue        
        
       end
       
        