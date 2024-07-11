function [P,A_Ac] = Inletfunction(u_Mach,u_alpha,inletdata,u_inlet)
%  INLETFUNCTION Brief summary of this function.
% 
% Detailed explanation of this function.
unique_alpha = unique(inletdata.Alpha)';
Inlet1index=find(inletdata.Inlet==1);
Inlet2index=find(inletdata.Inlet==2);
if u_inlet==1
    for counter = 1:5
        if u_alpha<=unique_alpha(counter)
            alpha_u = unique_alpha(counter);
            Ma_1=inletdata.M(find(inletdata.Alpha(Inlet1index)==unique_alpha(counter)));
            Po_1=inletdata.Po2bPoA(find(inletdata.Alpha(Inlet1index)==unique_alpha(counter)));
            Ac_1=inletdata.AabAc(find(inletdata.Alpha(Inlet1index)==unique_alpha(counter)));
            Pu=interp1(Ma_1,Po_1,u_Mach);
            Acu=interp1(Ma_1,Ac_1,u_Mach);
            break
        else
            continue
        end
    end
    for counter = 5:-1:1
        if u_alpha>=unique_alpha(counter)
            alpha_l = unique_alpha(counter);
            Ma_1=inletdata.M(find(inletdata.Alpha(Inlet1index)==unique_alpha(counter)));
            Po_1=inletdata.Po2bPoA(find(inletdata.Alpha(Inlet1index)==unique_alpha(counter)));
            Ac_1=inletdata.AabAc(find(inletdata.Alpha(Inlet1index)==unique_alpha(counter)));
            Pl=interp1(Ma_1,Po_1,u_Mach);
            Acl=interp1(Ma_1,Ac_1,u_Mach);
            break
        else
            continue
        end
    end
else
    for counter = 1:5
        if u_alpha<=unique_alpha(counter)
            alpha_u = unique_alpha(counter);
            Ma_2=inletdata.M(Inlet2index(inletdata.Alpha(Inlet2index)==unique_alpha(counter)));
            Po_2=inletdata.Po2bPoA(Inlet2index(inletdata.Alpha(Inlet2index)==unique_alpha(counter)));
            Ac_2=inletdata.AabAc(Inlet2index(inletdata.Alpha(Inlet2index)==unique_alpha(counter)));
            Pu=interp1(Ma_2,Po_2,u_Mach);
            Acu=interp1(Ma_2,Ac_2,u_Mach);
            break
        else
            continue
        end
    end
    for counter = 5:-1:1
        if u_alpha>=unique_alpha(counter)
            alpha_l = unique_alpha(counter);
            Ma_2=inletdata.M(Inlet2index(inletdata.Alpha(Inlet2index)==unique_alpha(counter)));
            Po_2=inletdata.Po2bPoA(Inlet2index(inletdata.Alpha(Inlet2index)==unique_alpha(counter)));
            Ac_2=inletdata.AabAc(Inlet2index(inletdata.Alpha(Inlet2index)==unique_alpha(counter)));
            Pl=interp1(Ma_2,Po_2,u_Mach);
            Acl=interp1(Ma_2,Ac_2,u_Mach);
            break
        else
            continue
        end
    end
end
P=interp1([alpha_l alpha_u],[Pl Pu], u_alpha);
A_Ac=interp1([alpha_l alpha_u],[Acl Acu], u_alpha);
end