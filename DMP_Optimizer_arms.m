function [func_value,c]=DMP_Optimizer7dot(q,E,q0i)
% OUTPUT:
% func_value: bình phương sai lech goc voi qoi
% c: sai lech vi tri
% INPUT:
% q: gia tri goc (rad)
% E: toa do diem cuoi (mm)
% q0i: gia tri goc cua diemr truoc
global canhTay cangTay rE
Omega=repmat([1 1 1 1 1 1 1],1,1);
 q1=q(:,1);q2=q(:,2);q3=q(:,3);q4=q(:,4);q5=q(:,5);q6=q(:,6);q7=q(:,7);
%  
% rE= [(607260097972743345*cos(q1)*cos(q3)*sin(q2)*sin(q4))/2251799813685248 - (2739417775299264423*cos(q1)*cos(q2))/9007199254740992 - (950399157810204315*cos(q4)*sin(q2))/72057594037927936 - (607260097972743345*sin(q1)*sin(q3)*sin(q4))/2251799813685248 - (607260097972743345*cos(q1)*cos(q2)*cos(q4))/2251799813685248 - (950399157810204315*cos(q2)*cos(q3)*sin(q4))/72057594037927936 - (4287356200788255021*sin(q2))/288230376151711744 - 1767038434150824319/288230376151711744;
%                                                                                                                                                                                                                                                                                                                                                             270*cos(q1)*sin(q3)*sin(q4) - 270*cos(q2)*cos(q4)*sin(q1) - (609*cos(q2)*sin(q1))/2 + 270*cos(q3)*sin(q1)*sin(q2)*sin(q4);
%   (2739417775299264423*sin(q2))/9007199254740992 - (4287356200788255021*cos(q1)*cos(q2))/288230376151711744 + (607260097972743345*cos(q4)*sin(q2))/2251799813685248 - (950399157810204315*sin(q1)*sin(q3)*sin(q4))/72057594037927936 - (950399157810204315*cos(q1)*cos(q2)*cos(q4))/72057594037927936 + (607260097972743345*cos(q2)*cos(q3)*sin(q4))/2251799813685248 + (950399157810204315*cos(q1)*cos(q3)*sin(q2)*sin(q4))/72057594037927936 + 1129053959934507997/9007199254740992]';
% 
% rE=[ (549*cos(q1)*cos(q3)*sin(q2)*sin(q4))/2 - (549*sin(q1)*sin(q3)*sin(q4))/2 - (549*cos(q1)*cos(q2)*cos(q4))/2 - (609*cos(q1)*cos(q2))/2;
%  (549*cos(q1)*sin(q3)*sin(q4))/2 - (549*cos(q2)*cos(q4)*sin(q1))/2 - (609*cos(q2)*sin(q1))/2 + (549*cos(q3)*sin(q1)*sin(q2)*sin(q4))/2;
%                                                    (609*sin(q2))/2 + (549*cos(q4)*sin(q2))/2 + (549*cos(q2)*cos(q3)*sin(q4))/2 + 251/2]';
% rE=[7*sin(q1) - (5733*cos(q1)*cos(q2))/20 - (193*cos(q1)*sin(q2))/20 - 266*sin(q1)*sin(q3)*sin(q4) - 266*cos(q1)*cos(q2)*cos(q4) - 266*cos(q1)*cos(q3)*sin(q2)*sin(q4)
%  266*cos(q1)*sin(q3)*sin(q4) - (5733*cos(q2)*sin(q1))/20 - (193*sin(q1)*sin(q2))/20 - 266*cos(q2)*cos(q4)*sin(q1) - 7*cos(q1) - 266*cos(q3)*sin(q1)*sin(q2)*sin(q4)
%                                                                      (193*cos(q2))/20 - (5733*sin(q2))/20 - 266*cos(q4)*sin(q2) + 266*cos(q2)*cos(q3)*sin(q4) + 143]';
%%Tinh toa do diem cuoi tu cong thuc rE tong quat
rEP=toa_do_diem_cuoi(q1,q2,q3,q4,q5);

xx=rEP(:,1)-E(:,1);
yy=rEP(:,2)-E(:,2);
zz=rEP(:,3)-E(:,3);
c=[xx yy zz];
%Objective Function

func_value=((q-q0i).^2).*Omega;
func_value=sqrt(sum(func_value,2));

end


