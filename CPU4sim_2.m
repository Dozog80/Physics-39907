function CPU4sim(prog)
%% Clear figure
clf;

%% Program
if(~exist('prog','var') || isempty(prog))
  
  [prog{1:16}]=deal('');
  
  % code (Division of A by B, result in 12, remainder in 13)

  %%% Subtract A and B (A - B), increment result -> store remainder -> repeat until negative. 

  prog{1}='LDA 14';  
  prog{2}='LDB 15'; 
  prog{3}='LDC 11'; % Store a 1 in C, used to increment result (in D)
  prog{4}='STA 13'; % Store remainders in 13
  prog{5}='SUB A B';  % Stored in A
  prog{6}='JPN 9';   %% Halt (since subtraction is negative, division is done)
  prog{7}='ADD D C';   % Increment result in D by 1 (This is stored in D!!!)
  prog{8}='STD 12';   % Store result in 12
  prog{9}='JMP 3'; % Loop over the subtractions
  prog{10}='HLT';          

  % data
  prog{12} = 1; % 1 stored in ADDR 11 to increment result each successful subtraction
  prog{15}= 9;          
  prog{16}= 11; % ADDR 15            
   
end

%% RAM
mem=defPanel('RAM','',[.71 .05 .288 .95]);

NRam=16;sep=.02;sz=(1-(NRam+1)*sep)/NRam;
wram=.95;
RAM(NRam)=struct('val',0,'name','');
for nr=NRam:-1:1
  RAM(nr).val=assembleASM(prog{nr});
  RAM(nr).name=sprintf('%02d',nr-1);
  pos=[1-sep-wram 1-(sep+sz)*nr wram sz];
  ramh(nr)=defREG(RAM(nr).name,RAM(nr).val,mem,pos);
  set(ramh(nr),'userdata','decode');
end

%% CPU
cpu=defPanel('CPU','',[.002 .5 .7 .5]);

%% Registers
reg=defPanel('Registers',cpu,[.01 .83 .98 .15],'R');

NR=4;sep=.005;sz=(1-(NR+1)*sep)/NR;
hr=.9;
R(NR)=struct('val',0,'name','');
for nr=NR:-1:1
  R(nr).val=randi(256)-1;
  R(nr).name='A'+nr-1;
  regh(nr)=defREG(R(nr).name,R(nr).val,reg,...
    [sep+(nr-1)*(sz+sep) (1-hr) sz hr]);
end

%% Control Unit
cu=defPanel('Control Unit',cpu,[.525 .02 .425 .8],'CU');

sep=.03;
sz=1/4;
% Instruction Register
ir=defPanel('Instruction Register',cu,[.1 1-sz .8 sz-sep],'IR');
IR=struct('val',randi(256)-1,'name','IR');
irh=defREG(IR.name,IR.val,ir,[.1 .1 .8 .8]);

% Instruction Decode
id=defPanel('Decoded Instruction',cu,[.1 1-2*sz .8 sz-sep],'IP');
NI=16;
INSL=setINST(NI);
idh=defREG('',0,id,[.1 .1 .8 .8]);
ID=struct('val',IR.val,'asm',decodeASM(IR.val));
updateID(ID,idh);

% Instruction Pointer
ip=defPanel('Instruction Pointer',cu,[.1 1-3*sz .8 sz-sep],'IP');
IP=struct('val',randi(16)-1,'name','IP');
iph=defREG(IP.name,IP.val,ip,[.1 .1 .8 .8],4);

% Flags
flags=defPanel('Flags',cu,[.1 1-4*sz .8 sz-sep],'F');
FLN={'O','N','Z'};
NF=length(FLN);
sep=.1;sz=.13;
for nf=NF:-1:1
  iFLN.(FLN{nf})=nf;
  FLAGS(nf)=struct('val',rand>.5,'name',FLN(nf));
  flagsh(nf)=defREG('',0,flags,[sep+sz*(nf-1) sep sz 1-2*sep]);
  c=[1 0 0]*~FLAGS(nf).val+[0 .8 0]*FLAGS(nf).val;
  set(flagsh(nf),'String',FLN(nf),'background',c,'horiz','center');
end

t={'T','F'};
for nf=2:-1:1
  lab(nf)=defREG('',0,flags,[sep+sz*(nf-1+NF+1) sep sz 1-2*sep]);
  c=[1 0 0]*(nf==2)+[0 .8 0]*(nf==1);
  set(lab(nf),'String',t(nf),'background',c,'horiz','center');
end


%% ALU
alu=defPanel('ALU',cpu,[.05 .22 .425 .6]);

% OPA
a=defPanel('Operand A',alu,[.05 1-.33 .65 .3],'OPA');
A=struct('val',randi(256)-1,'name','A');
ah=defREG(A.name,A.val,a,[.05 .1 .9 .8]);

% OPB
b=defPanel('Operand B',alu,[.05 1-.67 .65 .3],'OPB');
B=struct('val',randi(256)-1,'name','B');
bh=defREG(B.name,B.val,b,[.05 .1 .9 .8]);

% OP
op=defPanel('Op',alu,[.75 1-.55 .2 .4],'OP');
OPL={'+','-'};
OP=OPL{randi(2)};
oph=defREG('',0,op);
set(oph,'style','text','fontsize',15);
set(oph,'position',[.05 .1 .8 .8]);
set(oph,'String',OP,'horiz','center');

% OUT
out=defPanel('Output',alu,[.05 1-.99 .65 .3],'O');
val=randi(256)-1;
OUT=struct('val',val,'name','O');
outh=defREG(OUT.name,OUT.val,out,[.05 .1 .9 .8]);

% Over/under flow type oft
oft=defPanel('O/U flow',alu,[.72 1-.99 .27 .4],'OFT');
set(oft,'titlepos','centertop');
[ofth,obh]=defSwitch(2,oft,[.02 .02 .96 .96]);
set(ofth,'bordertype','none');
set(obh(1),'string','clip');
set(obh(2),'string','mod');

%% Fetch CPU

fet=uicontrol('style','pushbutton','units','normalized');
set(fet,'parent',cpu,'position',[.02 .01 .1 .1]);
set(fet,'string','Fetch');
set(fet,'callback',@fetch);

%% Decode CPU

dx=uicontrol('style','pushbutton','units','normalized');
set(dx,'parent',cpu,'position',[.02+.12 .01 .1 .1]);
set(dx,'string','Decode');
set(dx,'callback',@decode);

%% Execute CPU

exe=uicontrol('style','pushbutton','units','normalized');
set(exe,'parent',cpu,'position',[.02+.24 .01 .1 .1]);
set(exe,'string','Execute');
set(exe,'callback',@execute);

%% Reset CPU

rst=uicontrol('style','pushbutton','units','normalized');
set(rst,'parent',cpu,'position',[.02+.4 .01 .1 .1]);
set(rst,'string','Reset');
set(rst,'UserData',0);
set(rst,'callback',@reset);

%% Fetch/Decode CPU

rst=uicontrol('style','pushbutton','units','normalized');
set(rst,'parent',cpu,'position',[.07 .115 .1 .1]);
set(rst,'string','Fet/Dec');
set(rst,'callback',@fetdec);

%% Fetch/Decode/Excute CPU

rst=uicontrol('style','pushbutton','units','normalized');
set(rst,'parent',cpu,'position',[.07+.12 .115 .1 .1]);
set(rst,'string','Step');
set(rst,'callback',@step);

%% Fetch/Decode/Excute CPU

rst=uicontrol('style','pushbutton','units','normalized');
set(rst,'parent',cpu,'position',[.10+.24 .115 .1 .1]);
set(rst,'string','Run');
set(rst,'callback',@run);

%% Edit register callback
  function regedit(h,~)
    p=get(h,'parent');
    ptag=get(p,'tag');
    tag=get(h,'tag');
    switch ptag
      case 'RAM'
        in=get(h,'string');
        num=str2double(in);
        ad=str2double(tag);
        if(~isnan(num)) % if number
          RAM(ad+1).val=clip(num,0,255);
        else
          if(length(in)>2) % if 3 or longer
            inst=assembleASM(in);
            if(inst>=0) % if valid instruction
              RAM(ad+1).val=inst;
            end
          end
        end  
        updateREG(RAM(ad+1),ramh(ad+1));
      case 'R'
        in=get(h,'string');
        num=str2double(in);
        ad=tag-'A';
        if(~isnan(num))
          R(ad+1).val=clip(num,0,255);
        end
        updateREG(R(ad+1),regh(ad+1));
      otherwise
        updateCU();
        updateALU();
    end
  end

%% Run callback
  function err=run(h,e)
    err='';
    p=.1;
    set(rst,'userdata',0);
    while(~strcmp(err,'HLT') && get(rst,'userdata')==0)
      
      set(fet,'value',1);
      fetch(h,e);
      pause(p);
      set(fet,'value',0);
      
      set(dx,'value',1);
      decode(h,e);
      pause(p);
      set(dx,'value',0);
      
      set(exe,'value',1);
      err=execute(h,e);
      pause(p);
      set(exe,'value',0);
      
    end
    set(rst,'userdata',0);
  end

%% Fetch/Decode/Execute (step) callback
  function err=step(h,e)
    fetdec(h,e);
    err=execute(h,e);
  end

%% Fetch/Decode callback
  function err=fetdec(h,e)
    fetch(h,e);
    err=decode(h,e);
  end

%% Fetch callback
  function err=fetch(~,~)
    IR.val=RAM(IP.val+1).val;
    updateREG(IR,irh);
    IP.val=clip(IP.val+1,0,15);
    err=updateREG(IP,iph,4);
  end

%% decode callback
  function err=decode(~,~)
    ID.val=IR.val;
    ID.asm=decodeASM(ID.val);
    err=updateID(ID,idh);
  end

%% Execute callback
  function err=execute(~,~)
    err=executeIR(IR.val);
  end

%% Reset callback
  function err=reset(~,~)
    set(rst,'userdata',1);
    [R.val]=deal(0);
    A.val=0;
    B.val=0;
    OUT.val=0;
    OP='+';
    IR.val=15;
    IP.val=0;
    ID.val=15*16;
    ID.asm=decodeASM(ID.val);
    [FLAGS.val]=deal(false);
    updateREG(RAM,ramh);
    err=updateCPU;
  end

%% Update CPU
  function err=updateCPU()
    updateREG(R,regh);
    updateCU();
    err=updateALU();
  end

%% Update CU
  function err=updateCU()
    updateREG(IR,irh);
    updateREG(IP,iph,4);
    updateID(ID,idh);
    err=updateFLAGS(FLAGS,flagsh);
  end

%% Update ALU
  function err=updateALU()
    updateREG(A,ah);
    updateREG(B,bh);
    updateREG(OUT,outh);
    err=updateOP(OP,oph);
  end

%% Update REGISTER
  function err=updateREG(reg,h,bits)
    if(~exist('bits','var') || isempty(bits))
      bits=8;
    end
    N=length(reg);
    for n=1:N
      v=reg(n).val;
      if(v>2^bits-1);v=2^bits-1; end
      if(v<0);v=0; end
      if(strcmpi(get(h(n),'userdata'),'decode'))
        fmt=sprintf('%%s: %%s(%%0%dd) %%0%dX %%s',ceil(log10(2^bits-1)),ceil(log2(2^bits-1)/4));
        set(h(n),'String',sprintf(fmt,reg(n).name,dec2bin(v,bits),v,v,decodeASM(v)));
      else
        fmt=sprintf('%%s: %%s(%%0%dd)',ceil(log10(2^bits-1)));
        set(h(n),'String',sprintf(fmt,reg(n).name,dec2bin(v,bits),v));
      end
    end
    err=0;
  end

%% Update Flags
  function err=updateFLAGS(flg,h)
    N=length(flg);
    for n=1:N
      c=[1 0 0]*~flg(n).val+[0 .8 0]*flg(n).val;
      set(h(n),'background',c);
    end
    err=0;
  end

%% Update Decoded Instruction
  function err=updateID(reg,h)
    set(h,'string',reg.asm);
    err=0;
  end

%% Update ALU Operation
  function err=updateOP(oper,h)
    set(h,'String',oper);
    err=0;
  end

%% Decode Instruction (disassemble)
  function [asm,opc,ad]=decodeASM(in)
    opc=fix(in/16);
    ad=rem(in,16);
    nme=INSL(opc+1);
    asm=nme{:};
    if(opc>=1 && opc<=12)
      asm=[asm sprintf(' %02d',ad)];
    end
    if(opc>=13 && opc<=14)
      ad1=fix(ad/4)+'A';
      ad2=rem(ad,4)+'A';
      asm=[asm sprintf(' %c %c',ad1,ad2)];
    end
  end

%% Excute Instruction
  function [nme,opc,ad,err]=executeIR(in)
    opc=fix(in/16);
    ad=rem(in,16);
    nme=INSL{opc+1};
    if(opc>=1 && opc<=12)
      switch nme(1:2)
        case 'LD'                         % LOAD
          n=nme(3)-'A'+1;
          R(n).val=RAM(ad+1).val;
          updateREG(R(n),regh(n));
        case 'ST'                         % STORE
          n=nme(3)-'A'+1;
          RAM(ad+1).val=R(n).val;
          updateREG(RAM(ad+1),ramh(ad+1));
      end
      if(nme(1)=='J')                     % JUMP
        cd=nme(3);
        if((cd=='P')||...
            ((cd=='O') && FLAGS(iFLN.('O')).val)||...
            ((cd=='N') && FLAGS(iFLN.('N')).val)||...
            ((cd=='Z') && FLAGS(iFLN.('Z')).val))
          IP.val=ad;
          updateREG(IP,iph,4);
        end
      end
    end
    
    if(opc>=13 && opc<=14)                % ADD/SUB
      ad1=fix(ad/4);
      ad2=rem(ad,4);
      A.val=R(ad1+1).val;
      B.val=R(ad2+1).val;
      OP=OPL{opc-12};
      excuteALU();
      R(ad1+1).val=OUT.val;
      updateREG(R(ad1+1),regh(ad1+1));
    end
    if(opc==15)
      IP.val=clip(IP.val-1,0,15);
      err=updateREG(IP,iph,4);
    end
  end

%% Execute ALU
  function res=excuteALU()
    switch OP
      case '+'
        res=A.val+B.val;
      case '-'
        res=A.val-B.val;
    end
    
    [FLAGS.val]=deal(false);
    if(res>255);FLAGS(iFLN.('O')).val=true;end
    if(res<0);  FLAGS(iFLN.('N')).val=true;end
    if(res==0); FLAGS(iFLN.('Z')).val=true;end
    res=wrap(res,256,ofth);
    OUT.val=res;
    updateALU();
    updateFLAGS(FLAGS,flagsh);
  end

%% setup Register

  function h=defREG(name,val,parent,position,bits)
    % defReg <Default setup for register.>
    % Usage:: h=defReg(name,val,parent[0],positions[0])
    %
    
    % revision history:
    % 09/06/2023 Mark D. Shattuck <mds> defReg.m
    
    %% Create ui
    h=uicontrol('Style','edit');
    
    %% Parse Input
    if(~exist('parent','var') || isempty(parent))
    else
      set(h,'parent',parent);
    end
    
    set(h,'units','normalized');
    if(~exist('position','var') || isempty(position))
    else
      set(h,'Position',position);
    end
    
    if(~exist('bits','var') || isempty(bits))
      bits=8;
    end
    
    %% Main
    
    set(h,'fontname','monospaced');
    set(h,'HorizontalAlignment','left');
    
    MX=2^bits-1;
    val=clip(val,0,MX);
    fmt=sprintf('%%s: %%s(%%0%dd)',ceil(log10(MX)));
    set(h,'String',sprintf(fmt,name,dec2bin(val,bits),val));
    set(h,'callback',@regedit);
    tagname=sprintf('%s',name);
    set(h,'Tag',tagname);
  end
end

%% setup Panel

function h=defPanel(name,parent,position,tag)
% defPanel <Default setup for panel.>
% Usage:: h=defPanel(name,parent[0],positions[0],tag[name])
%

% revision history:
% 09/06/2023 Mark D. Shattuck <mds> defReg.m

%% Create ui
h=uipanel('Title',name);

%% Parse Input
if(~exist('parent','var') || isempty(parent))
else
  set(h,'parent',parent);
end

set(h,'units','normalized');
if(~exist('position','var') || isempty(position))
else
  set(h,'Position',position);
end

if(~exist('tag','var') || isempty(tag))
  tag=name;
end


%% Main

set(h,'Tag',tag);

end

%% setup Switch

function [h,bh]=defSwitch(num,parent,position)
% defSwitch <Default setup for switch panel.>
% Usage:: h=defSwitch(name,val,parent[0],positions[0])
%

% revision history:
% 09/05/2024 Mark D. Shattuck <mds> defSwitch.m

%% Create button group
h = uibuttongroup();

%% Parse Input
if(~exist('parent','var') || isempty(parent))
else
  set(h,'parent',parent);
end

set(h,'units','normalized');
if(~exist('position','var') || isempty(position))
else
  set(h,'Position',position);
end

%% Main
bh=zeros(1,num);
for n=1:num
  bh(n)= uicontrol('Style','Radio');
  set(bh(n),'parent',h);
  str=sprintf('opt %d',n);
  set(bh(n),'String',str);
  set(bh(n),'units','normalized');
  set(bh(n),'position',[.01 (num-n)/num .99 1/num]);
end
%set(h,'SelectionChangeFcn',@selcbk);
set(h,'SelectedObject',bh(1));
end

%% button change
% function selcbk(source,eventdata)
% disp(source);
% disp([eventdata.EventName,'  ',... 
%      get(eventdata.OldValue,'String'),'  ', ...
%      get(eventdata.NewValue,'String')]);
% disp(get(get(source,'SelectedObject'),'String'));
% end

%% Assemble ASM

function [inst,opc,ad,ad1,ad2]=assembleASM(asm)
ad=0;ad1=0;ad2=0;
iFLN=struct('O',1,'N',2,'Z',3);
if(isempty(asm))
  inst=0;
  opc=0;
  return;
end
if(isnumeric(asm))
  inst=clip(asm,0,255);
  opc=-1;
else
  nme=upper(asm(1:3));
  inst=true;
  switch nme(1:2)
    case 'LD'
      opc=2+nme(3)-'A'-1;
    case 'ST'
      opc=6+nme(3)-'A'-1;
    case 'JM'
      opc=10-1;
    case 'JP'
      opc=10-1+iFLN.(nme(3));
    otherwise
      inst=false;
  end
  if(inst)  % is LD/ST/J
    ad=str2double(asm(5:end));
    inst=opc*16+ad;
  else
    inst=true;
    switch nme
      case 'NOP'
        opc=0;
      case 'HLT'
        opc=15;
      otherwise
        inst=false;
    end
    if(inst) % is NOP/HLT
      inst=16*opc;
    else
      inst=true;
      switch nme
        case 'ADD'
          opc=14-1;
        case 'SUB'
          opc=15-1;
        otherwise
          inst=false;
      end
      if(inst) % is ADD/SUB
        ad1=asm(5)-'A';
        ad2=asm(7)-'A';
        inst=opc*16+ad1*4+ad2;
      else % unrecognized
        opc=-2;
        inst=15;
      end
    end
  end
end
end

%% Setup Instructions

function inst=setINST(N)
[inst{1:N}]=deal('NOP');
inst{N}='HLT';
%inst(1)='NOP'; %NOP No operation

inst{2}='LDA'; %LDA ADDR; Load RAM(ADDR) into reg A
inst{3}='LDB'; %LDB ADDR; Load RAM(ADDR) into reg B
inst{4}='LDC'; %LDC ADDR; Load RAM(ADDR) into reg C
inst{5}='LDD'; %LDD ADDR; Load RAM(ADDR) into reg D

inst{6}='STA'; %STA ADDR; Store reg A into RAM(ADDR)
inst{7}='STB'; %STA ADDR; Store reg B into RAM(ADDR)
inst{8}='STC'; %STA ADDR; Store reg C into RAM(ADDR)
inst{9}='STD'; %STA ADDR; Store reg D into RAM(ADDR)

inst{10}='JMP'; %JMP ADDR; Set Inst. poiter (IP) to ADDR
inst{11}='JPO'; %JZ  ADDR; Set IP to ADDR if Zero
inst{12}='JPN'; %JO  ADDR; Set IP to ADDR if Overflow
inst{13}='JPZ'; %JN  ADDR; Set IP to ADDR if Negative

inst{14}='ADD'; %ADD R1 R2; R1=R1+R2
inst{15}='SUB'; %SUB R1 R2; R1=R1-R2

%inst(16)='HLT'; %HLT; Decrement IP and stop
end

%% Clip (helper)

function out=wrap(in,max,h)
% Wraps in to N-bits based on current O/U flow model
% Usage:: out=wrapNbit(in,N)
%

% revision history:
% 09/06/2024 Mark D. Shattuck <mds>   wrapNbit.m

switch get(get(h,'SelectedObject'),'String')
  case 'clip'
    out=clip(in,0,max-1);
  case 'mod'
    out=mod(in,max);
end


end

%% Clip (helper)

function out=clip(in,lo,hi,ulo,uhi)
% clip Clips the values of in between lo and hi.
% Usage:: out=clip(in,lo,hi)
%
% Clips the values of in between lo and hi.

% revision history:
% 05/02/95 Mark D. Shattuck <mds>   clip.m

if(~exist('ulo','var')||isempty(ulo))
  ulo=lo;
end
if(~exist('uhi','var')||isempty(uhi))
  uhi=hi;
end

out=in;
out((out<lo))=ulo;
out((out>hi))=uhi;
end