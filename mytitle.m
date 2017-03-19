%ver 5.1
%19.03.2017
% Internal function for scattering patterns

function [ output ] = mytitle( ED,MD,EQ,MQ,EO,TD )

%%
switch ED
    case 1
    EDstr='Electric dipole, ';
    if MD+EQ+MQ+EO+TD==0
        EDstr='Electric dipole';
    end
        case 0
    EDstr='';
end

%%
switch MD
    case 1
        MDstr='Magnetic dipole';
        if EQ+MQ+EO+TD>0
            MDstr='Magnetic dipole, ';
        end
        if ED==1
            MDstr='Electric and magnetic dipoles';
            EDstr='';
        end
        if ED==1 && EQ+MQ+EO+TD>0
            MDstr='Electric and magnetic dipoles, ';
            EDstr='';
        end
    case 0
        MDstr='';
end

%%
switch TD
    case 1
        TDstr='Toroidal dipole';
        if EQ+MQ+EO>0
            TDstr='Toroidal dipole, ';
        end
        if ED==1 && MD~=1 && EQ+MQ+EO==0
            TDstr='Electric and toroidal dipoles';
            EDstr='';
        end
        if ED~=1 && MD==1 && EQ+MQ+EO==0
            TDstr='Magnetic and toroidal dipoles';
            MDstr='';
        end
        if ED==1 && MD==1 && EQ+MQ+EO==0
            TDstr='Electric, magnetic and toroidal dipoles';
            EDstr=''; MDstr='';
        end
        
        if ED==1 && MD~=1 && EQ+MQ+EO>0
            TDstr='Electric and toroidal dipoles, ';
            EDstr='';
        end
        if ED~=1 && MD==1 && EQ+MQ+EO>0
            TDstr='Magnetic and toroidal dipoles, ';
            MDstr='';
        end
        if ED==1 && MD==1 && EQ+MQ+EO>0
            TDstr='Electric, magnetic and toroidal dipoles, ';
            EDstr=''; MDstr='';
        end
    case 0
        TDstr='';
end
%%
switch EQ
    case 1
        EQstr='Electric quadrupole';
        if MQ+EO>0
            EQstr='Electric quadrupole, ';
        end
        if ED+MD+TD>0 
            EQstr='electric quadrupole';
        end
        if ED+MD+TD>0 && MQ+EO>0
            EQstr='electric quadrupole, ';
        end
    case 0
        EQstr='';
end

%%
switch MQ
    case 1
        MQstr='Magnetic quadrupole';
        if EO>0
            MQstr='Magnetic quadrupole, ';
        end
        if ED+MD+TD+EQ>0
            MQstr='magnetic quadrupole';
        end
        if ED+MD+TD+EQ>0 && EO>0 
            MQstr='magnetic quadrupole, ';
        end
        if EQ==1 && ED+MD+TD>0 && EO==0
            MQstr='electric and magnetic quadrupoles';
            EQstr='';
        end
        if EQ==1 && ED+MD+TD>0 && EO>0
            MQstr='electric and magnetic quadrupoles, ';
            EQstr='';                
        end
        if EQ==1 && ED+MD+TD==0 && EO==0
            MQstr='Electric and magnetic quadrupoles';
            EQstr='';                
        end
        if EQ==1 && ED+MD+TD==0 && EO>0
            MQstr='Electric and magnetic quadrupoles, ';
            EQstr='';                
        end
    case 0
        MQstr='';
end
            
%%
switch EO
    case 1
        EOstr='Electric octupole';
        if ED+MD+EQ+MQ+TD>0
            EOstr='electric octupole';
        end
        if ED+MD+TD+EQ+MQ==5
            EOstr='All multipoles';
            EDstr='';MDstr='';TDstr='';EQstr='';MQstr='';
        end
    case 0
        EOstr='';
end

output=[EDstr MDstr TDstr EQstr MQstr EOstr];


end


