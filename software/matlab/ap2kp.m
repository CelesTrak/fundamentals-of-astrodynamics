
% ------------------------------------------------------------------------------
%
%                           function ap2kp
%
%  this function converts ap to kp using cubic splines
%
%  author        : david vallado                  719-573-2600   4 aug  2005
%
%  revisions
%
%  inputs          description                    range / units
%    kpin        - kp 
%
%  outputs       :
%    apout       - ap
%
%  locals        :
%                - 
%
%  coupling      :
%
%  references    :
%    vallado       2004, 899-901
%
% kpout = ap2kp(apin);
% ------------------------------------------------------------------------------

function kpout = ap2kp(apin);
       
       ap = [-0.00001 -0.001 0 2 3 4 5 6 7 9 12 15 18 22 27 32 39 48 56 67 80 94 111 132 154 179 207 236 300 400 900];
       kp = [-0.66666667 -0.33333 0.0 0.33333 0.66667 1.0 1.33333 1.66667 2.0 2.33333 2.66667 ...
             3.0 3.33333 3.66667 4.0 4.33333 4.66667 5.0 5.33333 5.66667 ...  
             6.0 6.33333 6.66667 7.0 7.33333 7.66667 8.0 8.33333 8.66667  9.0 9.33333];

       kpout = 0.0;
       % now normalize the time to be 0, t1, t2, 1
%       t1 = 1.0/3.0;
%       t2 = 2.0/3.0;

       % find starting point in files based on input
       i = 1;
       while ((i<30) & (apin > ap(i)))
           i = i + 1;
         end;  
       i+1;
       
            if (i > 3)
               % -------- assign function points ---------
               p1 = kp(i - 2);
               p2 = kp(i - 1);
               p3 = kp(i );
               p4 = kp(i + 1);
               [ac0,ac1,ac2,ac3] = cubicspl(p1,p2,p3,p4);

               p1 = ap(i - 2);
               p2 = ap(i - 1);
               p3 = ap(i );
               p4 = ap(i + 1);
               [kac0,kac1,kac2,kac3] = cubicspl(p1,p2,p3,p4);
%p1
%p2
%p3
%p4
               % recover the original function values
               % use the normalized time first, but at an arbitrary interval
               [r1r,r1i,r2r,r2i,r3r,r3i] = cubic ( kac3,kac2,kac1,kac0-apin,'R' );
r(1) = r1r;
r(2)=r2r;
r(3)=r3r;
%               r = roots([kac3 kac2 kac1 kac0-apin])
%kac3
%kac2
%kac1
%kac0-apin
                if ((r(1) >= 0.00) & (r(1) <= 1.001))
                   apval = r(1);
                 else
                   if ((r(2) >= 0.00) & (r(2) <= 1.001))
                       apval = r(2);
                   else
                     if ((r(3) >= 0.00) & (r(3) <= 1.001))
                         apval = r(3);
                     else
                         apval = 0.0;  
                         fprintf(1,'error in root %11.7f %3i %11.7f %11.7f %11.7f \n',apin,i,r(1),r(2),r(3));
                     end;  
                   end; 
                 end;
%apin
%apval
                 kpout =  (ac3*apval^3  + ac2*apval^2  + ac1*apval + ac0);

%               fprintf(1,'ap2kp %11.7f  %11.7f \n',apval,kpout );
               
             end; % if i > 3
end                          

             
             
%  for i=0:27
%      ii = i*1.0/3.0;
%      kpin = ii;
%      apout = kp2ap(kpin);
%      kpout = ap2kp(apout);             
%      fprintf(1,'kp %11.7f  %11.7f  %11.7f \n',ii, apout, kpout);
%  end

