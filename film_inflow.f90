
module subprog
  implicit none
  !入力ぱらめーた
  !file_out_span
  integer::file_out_span=10
  !薄膜ぱらめーた
  integer,parameter::i_begin=24,i_end=25
  !25~26
  integer,parameter::j_begin=16,j_end=36
  real(8)::Voltage_film=-10.0d0
  !りある位置= (mesh-1)*dx

  
  
  !計算領域[m]
  real(8)::XX=50.0d-2,YY=50.0d-2

  !電子温度,電子,ion密度[3D] , 2D密度
  real(8)::Te=1.0d0             ![eV]
  real(8)::ne=(1.0d06)*1.0d06      ![m^-3] 
  real(8)::ne2d !=ne**(2.0/3.0) ![m^-2]
                !ne2d=ne??


  
  !drift速度[m/s]
  !real(8),parameter::v_drift=4.70d5

  
  !入力条件から一意的に決定されるパラメータ
  !ぷらずま角速度,でばい長
  real(8)::wpe             ![rad/s]
  real(8)::ld              ![m]
  !めっしゅ幅
  !CFL条件より
  !dx,dy<λd*10/3=2.46E-2[m]
  !dt<0.2/wpe=3.55E-9[ns]
  real(8)::dx=1.0d-2,dy=1.0d-2
  real(8)::dt=3.0d-09



  
  !!物理定数
  real(8)::EPS0=8.854187817d-12 ![F/m]誘電率
  real(8)::Q0=1.602176462d-19   ![C]素電荷量
  real(8)::kB=1.380649d-23      ![J/K]ぼるつまん
  real(8)::SM0=9.1091d-31       ![kg]電子の静止質量
  real(8)::SMXe=2.18017d-25   ![kg]Xe相対静止質量
  real(8)::pi=2*acos(0.0d0)     !!pi=Π



  !######################
  !######################
  !######################
  !!超粒子parameter
  real(8)::Q0_super,SM0_super  !超粒子一個あたりの電荷、質量
  real(8)::         SMXe_super !ion超粒子質量 電荷は同じ

  !1cellあたりの超粒子数
  !これを変更すると、超粒子質量が変わるので、
  !maxwellian生成るーちんの速度step幅 dv(local変数)も適当な値に変えてね
  integer,parameter::N_super_per_cell=16

  integer::N_per_cell          !1格子あたりの実粒子数
  integer::N_super             !超粒子一個あたりの実粒子数
  integer::N_super_number      !計算領域全超粒子数

  
  real(8)::TT
  integer::N !!=int(ne2d*XX*YY) [個] !!使うことはまず無い
  integer::N_count !!計算領域の総粒子数をcountする。


  !!空間３軸,粒子数,時間軸counter
  integer::i,j,k,l,t
  
  !空間めっしゅ数
  integer::NX,NY
  !!時間めっしゅの最大計算回数
  integer::NT
  !!経過時間step回数をかうんと
  integer::tstep_count




  !粒子情報、格子情報
  !↓↓↓↓↓↓↓↓↓↓↓↓
  !粒子情報
  !!それぞれ位置,その方向への速さ
  !一様分布させる粒子
  real(8),allocatable::location_x(:),location_y(:)
  real(8),allocatable::velocity_x(:),velocity_y(:)
  !電場
  real(8),allocatable::EX_particle(:),EY_particle(:)
  
  !                      領域内、外存在識別子 と 電子、いおん識別子
  !           .false.=境界条件含め、領域外へ飛び出た粒子
  !                                       .flase.=電子 
  logical,allocatable::lg_array_loc(:),lg_array_charge(:)

  !境界から流入させる粒子ver        (:,:)=(ΔNall*4+1,NT)
  real(8),allocatable::flow_loc_x(:,:),flow_loc_y(:,:)
  real(8),allocatable::flow_vel_x(:,:),flow_vel_y(:,:)
  real(8),allocatable::flow_EXP(:,:),flow_EYP(:,:)
  logical,allocatable::flow_lg_loc(:,:),flow_lg_cha(:,:)

  !Δt,1meshあたりに流入させる電子数
  integer::delta_Nsuper
  real(8)::delta_Nsuper_float
  !片側から流入させる粒子数
  !delta_Nall=delta_Nsuper*NX orNY
  integer::delta_Nall
                        !負方向から流入    !正方向から流入
  real(8),allocatable::vx_inflow_neg(:),vx_inflow_pos(:)
  real(8),allocatable::vy_inflow_neg(:),vy_inflow_pos(:)
  real(8),allocatable::v_inflow_thermal(:)


  

  !↑
  !重み関数などを駆使して、変換
  !↓
 
  !!格子情報
  real(8),allocatable::rho(:,:)
  real(8),allocatable::charge(:,:)
  real(8),allocatable::EX(:,:),EY(:,:),EX_half(:,:),EY_half(:,:)
  real(8),allocatable::phi(:,:)


  !MKSC単位 規格化定数
  real(8)::x_normal,m_normal,t_normal,Q_normal
  !組み立て単位
  real(8)::v_normal,charge_normal,rho_normal,phi_normal,E_normal,EPS_normal

  !maxwellian
  real(8)::vth_p3d

  
  




contains
  subroutine input_allocate
    real(8)::sheath,bohm_pot
    character(len=30)::fname_input_parameter
    
    !3dぴーく値熱速度
    vth_p3d=sqrt(2*Q0*Te/SM0)
    !計算mesh数
    NX=nint(XX/(dx))
    NY=nint(YY/(dy))

    !でばい[m]
    ld=sqrt(EPS0*kB*Te*11600/(ne*Q0**2))
    !角周波数[s]
    wpe=sqrt((ne*Q0**2)/(SM0*EPS0))

    !2D電子密度[1/m^2]
    ne2d=ne

    
    sheath=ld/sqrt(0.97d0)*(Q0*abs(Voltage_film)/(kB*Te*11600))**(0.75d0)
    bohm_pot=(kB*Te*11600)/(2*Q0)
    
    print'("#計算ぱらめーたを表示します")'
    print'("#     計算領域")'
    print'("#X方向[cm]  |   Y方向[cm]")'
    print'("#"F8.3,5X,F8.3)',XX*1.0d2,YY*1.0d2
    print'()'

 
    print'()'
    print'("# 電子密度[cm^-2]  |  電子温度[eV]")'
    print'("#"E10.3,10X,E10.3)',ne2d*1.0d-4,Te
    print'("# でばい長[cm]    |   角周波数[rad/s]")'
    print'("#",F10.3,10X,E15.5)',ld*1.0d2,wpe
    print'()'

    print'("# ion温度=0eVのしーす厚[cm]")'
    print'(E15.5)',sheath*1.0d2
    print'("# Bohm電位[V]")'
    print'(E15.5)',bohm_pot
    
    print'("#格子間隔,計算領域はCFL条件に従う")'
    print'("#dt<0.2/wpe  ,  dx,dy<3.3*ld")'
    print'("#dt=",F5.3,"[ns] | dx=",F5.3,"[cm] | dy=",F5.3"[cm]")'&
         ,dt*1.0d9,dx*1.0d2,dy*1.0d2
    print'()'
    


    !N_super_per_cell はgrobal変数にて定義
    !N_super_per_cell 1cellあたりの粒子数
    N_per_cell=ne2d*dx*dy       !計算格子1cellあたりの実粒子数
    
    N_super=int(ne2d*dx*dy/N_super_per_cell)   !超粒子1個あたりの実粒子数
                       
    N_super_number=N_super_per_cell*NX*NY
                               
    Q0_super=Q0*dble(N_super)
    SM0_super=SM0*dble(N_super)
    SMXe_super=SMXe*dble(N_super)
    print'("#超粒子parameter")'
    print'("#1meshあたりの超粒子数")'
    print'("#電子 + いおん =")'
    print'("#",I2,"+",I2)',N_super_per_cell/2,N_super_per_cell/2
    print'("#超粒子1個あたりの実粒子数",E15.5)',dble(N_super)
    print'("#超粒子電荷",E15.5)',Q0_super
    print'("#超粒子電子質量",E15.5)',SM0_super
    print'("#超粒子いおん質量",E15.5)',SMXe_super
    print'("#体系内の全超粒子数",I15.10)',N_super_number
    
    print'("###計測時間最大値を入力してください。[ns]###")'
    read*,TT
    TT=TT*1.0d-09     
    NT=nint(TT/(dt)) !!時間step回数
    

    allocate(location_x(N_super_number),location_y(N_super_number))
    allocate(velocity_x(N_super_number),velocity_y(N_super_number))
    allocate(EX_particle(N_super_number),EY_particle(N_super_number))
    allocate(lg_array_loc(N_super_number),lg_array_charge(N_super_number))
    
    
    allocate(rho(NX+2,NY+2),charge(NX+2,NY+2))
    allocate(EX(NX,NY),EY(NX,NY),EX_half(NX+2,NY+2),&
         EY_half(NX+2,NY+2),phi(NX+2,NY+2))



    
    write(fname_input_parameter,'(a)') 'imput_parameter.dat' 
    open(100,file=fname_input_parameter , status='replace')
    write(100,*) '##########  input parameter  ###############'
    write(100,*) 'mesh size[m],[s],'
    write(100,*) 'dx  dy  dt'
    write(100,'(2F10.3,E15.5)') dx,dy,dt

    write(100,*) '計算領域[m],計算時間[s]'
    write(100,*) 'X方向     Y方向    cal time'
    write(100,'(2F10.3,E15.5)') XX,YY,TT

    write(100,*) 'plasma palameter'
    write(100,*) '電子温度[eV],電子密度[1/m^3],質量[kg]'
    write(100,'(F10.3,2E15.5)') Te,ne,SM0
    
    write(100,*) 'ion温度[eV],ion密度[1/m^3],ion質量[kg]'
    write(100,'(a,2E15.5)') '0.0',ne,SMXe

    write(100,*) '角速度[rad/s],でばい長[m]'
    write(100,'(2E15.5)') wpe,ld

    write(100,*) '薄膜挿入位置[m]'
    write(100,*) 'x方向                        y方向'
    write(100,'(F10.5,a,F10.5,a,F10.5,a,F10.5)') i_begin*dx, '~' ,i_end*dx,&
         ',',j_begin*dy, '~' ,j_end*dy
    
    write(100,*) '薄膜電位[V]'
    write(100,'(F10.3)') Voltage_film
    
    write(100,*) 'ion温度0[eV]のしーす厚[m]'
    write(100,'(E10.3)') sheath

    write(100,*) ' Bohm電位[V]'
    write(100,'(E10.3)') bohm_pot
    
    write(100,*) 'Bohm速度[m/s]'
    write(100,'(E10.3)') sqrt(kB*Te*11600/SMXe)
    

    write(100,*) '超粒子parameter'
    write(100,*) '1cellあたりの超粒子数(1cellあたりの経算量)'
    write(100,'(I10.3)') N_super_per_cell
    write(100,*) '1超粒子あたりの実粒子数'
    write(100,'(I10.3)') N_super
    write(100,*) '総計算粒子数'
    write(100,'(I10.3)') N_super_number
    write(100,*) 
    close(100)


  end subroutine input_allocate


  subroutine passing_to_main(NT_sub,file_out_span_sub,&
       Te_sub,SM0_super_sub,SMXe_super_sub)
    !時間step数
    integer,intent(out)::NT_sub
    !file output span(ふぁいるのoutput間隔)
    integer,intent(out)::file_out_span_sub
    !電子温度
    real(8),intent(out)::Te_sub
    !それぞれ超粒子質量
    real(8),intent(out)::SM0_super_sub,SMXe_super_sub

    !引渡し変数
    NT_sub=NT
    file_out_span_sub=file_out_span
    Te_sub=Te
    SM0_super_sub=SM0_super
    SMXe_super_sub=SMXe_super
  end subroutine passing_to_main
  
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  
  
  subroutine charge_decided
    
    do l=1,N_super_number,2
       lg_array_charge(l)=.true. !正電荷
    end do

    do l=2,N_super_number,2
       lg_array_charge(l)=.false. !負電荷
    end do
    

  end subroutine charge_decided
  


  

  subroutine normalization
    character(len=30)::fname_normalize_para
    

    !規格化定数の定義
    x_normal=dx !(=dy)  ![m]
    m_normal=SM0_super  ![kg]
    t_normal=dt         ![s]
    Q_normal=Q0_super   ![C]

    v_normal=x_normal/t_normal
    charge_normal=Q_normal
    rho_normal=Q_normal/(x_normal**2)
    phi_normal=(m_normal*(x_normal**2))/(Q_normal*(t_normal**2))
    E_normal=(m_normal*x_normal)/(Q_normal*(t_normal**2))
    EPS_normal=((t_normal**2)*(Q_normal**2))/(m_normal*(x_normal**3))

    write(fname_normalize_para,'(a)') 'normalize_parameter.dat' 
    open(100,file=fname_normalize_para , status='replace')
    write(100,*) '##########  normalize parameter  ###############'
    write(100,*) '[m],[kg],[s],[C]'
    write(100,'(4E15.4)') x_normal,m_normal,t_normal,Q_normal

    write(100,*) '[m/s],[C/m^2],[V],[V/m],[F/m]'
    write(100,'(5E15.4)') v_normal,rho_normal,phi_normal,E_normal,EPS_normal
    
    close(100)


    
    !規格化
    ![m]
    dx=dx/x_normal
    dy=dy/x_normal
    
    XX=XX/x_normal
    YY=YY/x_normal
    
    location_x(:)=location_x(:)/x_normal
    location_y(:)=location_y(:)/x_normal

    
    
    ne=ne*x_normal**(-3.0d0)
    ne2d=ne2d*x_normal**(-2.0d0)
    

    ![kg]
    SM0_super=SM0_super/m_normal
    SMXe_super=SMXe_super/m_normal
    
    ![s]
    dt=dt/t_normal

    ![C]
    Q0_super=Q0_super/Q_normal
    charge(:,:)=charge(:,:)/charge_normal

    ![m/s]
    velocity_x(:)=velocity_x(:)/v_normal
    velocity_y(:)=velocity_y(:)/v_normal
    !流入時の速度
    vx_inflow_neg(:)=vx_inflow_neg(:)/v_normal
    vx_inflow_pos(:)=vx_inflow_pos(:)/v_normal
    vy_inflow_neg(:)=vy_inflow_neg(:)/v_normal
    vy_inflow_pos(:)=vy_inflow_pos(:)/v_normal
    v_inflow_thermal(:)=v_inflow_thermal(:)/v_normal


    ![F/m]=[s^2・C^2/kg・m^3]
    EPS0=EPS0/EPS_normal

    
    !以下はやらなくてもよいが、念のため
    ![C/m^2]
    rho(:,:)=rho(:,:)/rho_normal

    ![V]=[kg・m^2/s^2・C]
    phi(:,:)=phi(:,:)/phi_normal
    Voltage_film=Voltage_film/phi_normal

    ![V/m]=[kg・m/s^2・C]
    EX(:,:)=EX(:,:)/E_normal
    EY(:,:)=EY(:,:)/E_normal
    EX_particle(:)=EX_particle(:)/E_normal
    EY_particle(:)=EY_particle(:)/E_normal

    
  end subroutine normalization
  


















  



















  

  



  subroutine initial_conditions

    real(8)::x(N_super_number),y(N_super_number)
    real(8)::hanteix,hanteiy
    
    
    integer::counter
    integer::clock_time,seed_size
    real(8)::RN(N_super_number)
    integer,allocatable::seed(:)


    call system_clock(clock_time)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed=clock_time
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(RN(:))



    print'("Generating uniform distribution now")'

    do i=1,NX
       x(N_super_per_cell*(i-1)+1:N_super_per_cell*i)=dx*(RN(N_super_per_cell*(i-1)+&
            1:N_super_per_cell*i)+(i-1))
    end do
    
    call system_clock(clock_time)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed=clock_time
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(RN(:))

    do j=1,NY
       y(N_super_per_cell*(j-1)+1:N_super_per_cell*j)=dy*(RN(N_super_per_cell*(j-1)+&
            1:N_super_per_cell*j)+(j-1))
    end do

    counter=1
    do j=1,NY
       do i=1,NX
          location_x(N_super_per_cell*(counter-1)+1:N_super_per_cell*counter)&
               =x(N_super_per_cell*(i-1)+1:N_super_per_cell*i)
          location_y(N_super_per_cell*(counter-1)+1:N_super_per_cell*counter)&
               =y(N_super_per_cell*(j-1)+1:N_super_per_cell*j)
          
          counter=counter+1
       end do
    end do

    allocate(flow_loc_x(delta_Nall*4,NT),flow_loc_y(delta_Nall*4,NT),&
         flow_vel_x(delta_Nall*4,NT),flow_vel_y(delta_Nall*4,NT),&
         flow_EXP(delta_Nall*4,NT),flow_EYP(delta_Nall*4,NT),&
         flow_lg_loc(delta_Nall*4,NT),flow_lg_cha(delta_Nall*4,NT))

    
  end subroutine initial_conditions










  
!以降時間るーぷで囲まれた計算
  subroutine counter
    tstep_count=tstep_count+1
    print'("##経過時間",F15.1,"[ns]")',3*dble(tstep_count)
  end subroutine counter
  


  subroutine initial_conditions_inflow
    integer::counter
    integer::clock_time,seed_size
    real(8)::RN(delta_Nall)
    integer,allocatable::seed(:)

    real(8)::loc_x(delta_Nall),loc_y(delta_Nall)
    
    
    call system_clock(clock_time)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed=clock_time
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(RN(:))

    !流入粒子 一様分布routine
    !do i=1,NX
    !   loc_x(delta_Nsuper*(i-1)+1:delta_Nsuper*i)=(RN(delta_Nsuper*(i-1)+&
    !        1:delta_Nsuper*i)+(i-1))
    !end do
    !do j=1,NY
    !   loc_y(delta_Nsuper*(j-1)+1:delta_Nsuper*j)=(RN(delta_Nsuper*(j-1)+&
    !        1:delta_Nsuper*j)+(j-1))
    !end do

    
    
       
       
    

    !速度代入
    !-xから
    flow_vel_x(1:delta_Nall,tstep_count)=vx_inflow_neg(:)
    flow_vel_y(1:delta_Nall,tstep_count)=v_inflow_thermal(:)

    !+xから
    flow_vel_x(delta_Nall+1:2*delta_Nall,tstep_count)&
         =vx_inflow_pos(:)
    flow_vel_y(delta_Nall+1:2*delta_Nall,tstep_count)&
         =v_inflow_thermal(:)

    !-yから
    flow_vel_x(2*delta_Nall+1:3*delta_Nall,tstep_count)&
         =v_inflow_thermal(:) 
    flow_vel_y(2*delta_Nall+1:3*delta_Nall,tstep_count)&
         =vy_inflow_neg(:)
    
    !+yから
    flow_vel_x(3*delta_Nall+1:4*delta_Nall,tstep_count)&
         =v_inflow_thermal(:)
    flow_vel_y(3*delta_Nall+1:4*delta_Nall,tstep_count)&
         =vy_inflow_pos(:)






    
    !位置代入
    !-x
    call system_clock(clock_time)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed=clock_time
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(RN(:))
    do l=1,delta_Nall
       loc_y(l)=RN(l)*YY
    end do
    flow_loc_x(1:delta_Nall,tstep_count)=flow_vel_x(1:delta_Nall,tstep_count)*dt*RN(:)
    flow_loc_y(1:delta_Nall,tstep_count)=loc_y(1:delta_Nall)
    

    !+x
    call system_clock(clock_time)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed=clock_time+1
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(RN(:))
    do l=1,delta_Nall
       loc_y(l)=RN(l)*YY
    end do
    flow_loc_x(delta_Nall+1:2*delta_Nall,tstep_count)&
         =XX+(flow_vel_x(delta_Nall+1:2*delta_Nall,tstep_count)*dt*RN(:))
    flow_loc_y(delta_Nall+1:2*delta_Nall,tstep_count)=loc_y(1:delta_Nall)

    
    !-y
    call system_clock(clock_time)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed=clock_time+2
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(RN(:))
    do l=1,delta_Nall
       loc_x(l)=RN(l)*XX
    end do
    flow_loc_x(2*delta_Nall+1:3*delta_Nall,tstep_count)=loc_x(1:delta_Nall)
    flow_loc_y(2*delta_Nall+1:3*delta_Nall,tstep_count)&
         =flow_vel_y(2*delta_Nall+1:3*delta_Nall,tstep_count)*dt*RN(:)



    
    !+y
    call system_clock(clock_time)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed=clock_time+3
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(RN(:))
    do l=1,delta_Nall
       loc_x(l)=RN(l)*XX
    end do
    flow_loc_x(3*delta_Nall+1:4*delta_Nall,tstep_count)=loc_x(1:delta_Nall)
    flow_loc_y(3*delta_Nall+1:4*delta_Nall,tstep_count)&
         =YY+(flow_vel_y(3*delta_Nall+1:4*delta_Nall,tstep_count)*dt*RN(:))




    !print'("        vx       |       vy      |       x       |        y")'
    !do l=0*delta_Nall+1,1*delta_Nall
    !   print'(4F15.5)',flow_vel_x(l,tstep_count)*v_normal,&
    !        flow_vel_y(l,tstep_count)*v_normal,&
    !        flow_loc_x(l,tstep_count)*x_normal,flow_loc_y(l,tstep_count)*x_normal
    !end do
    
   
  end subroutine initial_conditions_inflow







  subroutine C_density_distribute
    !!めっしゅ区間の密度 
    !!密度格納配列
    real(8)::hanteix,hanteiy
    real(8)::S1,S2,S3,S4
    charge(:,:)=0
    rho(:,:)=0

    !初期一様分布粒子
    do i=1,NX+1
       do j=1,NY+1
          do l=1,N_super_number
             hanteix=location_x(l)!!粒子のx座標
             hanteiy=location_y(l)!!粒子のy座標
             
             if(dx*(i-1)<=hanteix.and.hanteix<=dx*(i-1)+dx&
                  .and.dy*(j-1)<=hanteiy.and.hanteiy<=dy*(j-1)+dy)then
                S4=(hanteix-dx*(i-1))*(hanteiy-dy*(j-1))
                S3=(dx*i-hanteix)*(hanteiy-dy*(j-1))
                S2=(hanteix-dx*(i-1))*(dy*j-hanteiy)
                S1=(dx*i-hanteix)*(dy*j-hanteiy)
                !!対角の面積比から、電荷量を分配[C]
                if(lg_array_charge(l))then                      
                   charge(i,j)=charge(i,j)+Q0_super*S1/dx/dy
                   charge(i+1,j)=charge(i+1,j)+Q0_super*S2/dx/dy
                   charge(i,j+1)=charge(i,j+1)+Q0_super*S3/dx/dy
                   charge(i+1,j+1)=charge(i+1,j+1)+Q0_super*S4/dx/dy
                else if(.not.lg_array_charge(l))then    
                   charge(i,j)=charge(i,j)-Q0_super*S1/dx/dy
                   charge(i+1,j)=charge(i+1,j)-Q0_super*S2/dx/dy
                   charge(i,j+1)=charge(i,j+1)-Q0_super*S3/dx/dy
                   charge(i+1,j+1)=charge(i+1,j+1)-Q0_super*S4/dx/dy
                end if                   
             end if
          end do
       end do
    end do


    !流入粒子
    do i=1,NX+1
       do j=1,NY+1
          do l=1,4*delta_Nall
             do t=1,tstep_count
                hanteix=flow_loc_x(l,t)
                hanteiy=flow_loc_y(l,t)
                if(dx*(i-1)<=hanteix.and.hanteix<=dx*(i-1)+dx&
                     .and.dy*(j-1)<=hanteiy.and.hanteiy<=dy*(j-1)+dy)then
                   S4=(hanteix-dx*(i-1))*(hanteiy-dy*(j-1))
                   S3=(dx*i-hanteix)*(hanteiy-dy*(j-1))
                   S2=(hanteix-dx*(i-1))*(dy*j-hanteiy)
                   S1=(dx*i-hanteix)*(dy*j-hanteiy)

                !正負電荷判定の必要があれば追加   
                   charge(i,j)=charge(i,j)-Q0_super*S1/dx/dy
                   charge(i+1,j)=charge(i+1,j)-Q0_super*S2/dx/dy
                   charge(i,j+1)=charge(i,j+1)-Q0_super*S3/dx/dy
                   charge(i+1,j+1)=charge(i+1,j+1)-Q0_super*S4/dx/dy
                end if
             end do   
          end do
       end do
    end do

    
    !!電荷密度計算[C/cm^-2]
    rho(:,:)=charge(:,:)/dx/dy
    
  end subroutine C_density_distribute


  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  subroutine poisson_Gauss_Seidel
    !******************
    !  ∇・∇φ=-ρ/ε0   *
    !******************
    !!dx=dyのもとで下記の式(差分化)が成り立つことに注意
    integer::loop
    real(8)::pre_phi,Current_phi,Max_phi
    real(8)::MAX_error,Current_error
    real(8)::wsor,lamda
    real(8)::pre(NX+2,NY+2)
    
    lamda=0.5d0*(cos(pi/NX)+cos(pi/NY))
    wsor=2/(1-sqrt(lamda**2))

    !   #image#
    !               neumann
    !        ###################
    !        #!!!!!!!!!!!!!!!!!#
    !        #!               !#
    !        #!  calculation  !#
    ! neumann#!       area    !#  neumann
    !        #!               !#
    !        #!!!!!!!!!!!!!!!!!#
    !        ###################
    !             neumann

    
    !!電位計算
    !!最大電位で規格化した電位誤差量が、一定値を下回った時ガウス＝ザイダルの反復法を終了
    loop=0 !!loop counter
    Max_phi=1.0d-10
    Current_phi=0.0d0

    pre(:,:)=0
    phi(:,:)=0

    loop0:do 
       Max_error=0.0d0
       Current_error=0.0d0
       loop1:do i=1,NX+1
          loop2:do j=1,NY+1
             !!ポアソン差分化
             pre_phi=phi(i,j)!!前回分の電位との差を計算するため、前回の値を一時的に保存

             !##のいまん型境界条件
             !!計算領域端では、電位勾配を0に保つ
             !糊代を付与したが、内側の値を参照するとなお良い
             if(j==1)then
                phi(i,j)=phi(i,j+1)
             end if
             
             if(j==NY+1)then
                phi(i,j)=phi(i,j-1)
             end if
             
             if(i==1)then
                phi(i,j)=phi(i+1,j)
             end if
             
             if(i==NX+1)then
                phi(i,j)=phi(i-1,j)
             end if



             !if(j==1)then
             !   phi(i,j)=phi(i,j+1)
             !else if(j==NY+1)then
             !   phi(i,j)=phi(i,j-1)
             !else if(i==1)then
             !   phi(i,j)=phi(i+1,j)
             !else if(i==NX+1)then
             !   phi(i,j)=phi(i-1,j)

                !物体付近の特殊solver
             !else if(i==i_begin-1.and.j_begin-1<=j.and.j<=j_end+1)then
             !   phi(i,j)=0.25d0*(2*phi(i+1,j)+phi(i,j-1)+phi(i,j+1)&
             !        +rho(i,j)*dx*dy/EPS0)

             !else if(i==i_end+1.and.j_begin-1<=j.and.j<=j_end+1)then
             !   phi(i,j)=0.25d0*(2*phi(i-1,j)+phi(i,j-1)+phi(i,j+1)&
             !        +rho(i,j)*dx*dy/EPS0)              
                
                
                !薄膜挿入位置(固定電位)
             if(i_begin<=i.and.i<=i_end.and.&
                  j_begin<=j.and.j<=j_end)then
                phi(i,j)=Voltage_film
             else
                !sor法 加速ぱらめーたをくわえた
                !phi(i,j)=pre(i,j)+0.25d0*wsor*(rho(i,j)*dx*dy/EPS0+&
                !     pre(i+1,j)+pre(i-1,j)+pre(i,j+1)+pre(i,j-1)&
                !     -4.0*pre(i,j))
                
                !がうすざいだる
                phi(i,j)=0.25d0*(rho(i,j)*dx*dy/EPS0+&
                     phi(i+1,j)+phi(i-1,j)+phi(i,j+1)+phi(i,j-1))

             end if
                   
             Current_phi=phi(i,j)
             Current_error=abs(pre_phi-Current_phi)!!前回との誤差量
             !!最大電位を更新
             if(Max_phi<Current_phi)then
                Max_phi=Current_phi
             end if
             !!最大誤差量を更新
             if(Max_error<Current_error)then
                MAX_error=Current_error
             end if
          end do loop2
       end do loop1
       
       loop=loop+1
       pre(:,:)=phi(:,:)
       
       if((Max_error/Max_phi)<=1.0d-10)exit loop0
       if(loop==1000000)exit loop0
                              !反復における,許容誤差量,許容回数
    end do loop0

    if(loop==1000000)then
       print'("収束しませんでした")'
    end if
    phi(NX+2,:)=phi(NX+1,:)
    phi(:,NY+2)=phi(:,NY+1)
    print'("##ガウス＝ザイデル反復回数:",I10)',loop
    
  end subroutine poisson_Gauss_Seidel


  subroutine phidis_output
    character(len=30)::fname_phi

    
    write(fname_phi,'(a,I7,a)') 'phit=',nint(3*dble(tstep_count)),'[ns].csv'
    
    open(100,file=fname_phi , status='replace')

    do i=1,NX+1
       do j=1,NY+1
          write(100,'(E15.3,",",E15.3,",",E20.10)')&
               dx*(i-1)*x_normal*1.0d2,&
               dy*(j-1)*x_normal*1.0d2,&
               phi(i,j)*phi_normal
       end do
       write(100,*) !xの値が変更する際、改行しないとgnuplotは正しく読み込んでくれない。
    end do
    close(100)
    
  end subroutine phidis_output
  
  !????????????????????????????????
  !????????????????????????????????

  subroutine E_field
    real(8)::hanteix,hanteiy
    real(8)::S1,S2,S3,S4
    
    EX(:,:)=0
    EY(:,:)=0
    EX_half(:,:)=0
    EY_half(:,:)=0
    EX_particle(:)=0
    EY_particle(:)=0
    flow_EXP(:,:)=0
    flow_EYP(:,:)=0
    

    do j=1,NY+1
       do i=1,NX+1
          if(i==1)then
             EX(i,j)=-(-3*phi(i,j)+4*phi(i+1,j)-phi(i+2,j))/(2*dx)
          else if(i==NX+1)then
             EX(i,j)=-(3*phi(i,j)-4*phi(i-1,j)+phi(i-2,j))/(2*dx)
          else
             EX(i,j)=-(phi(i+1,j)-phi(i-1,j))/(2*dx)
          end if
       end do
    end do


    
    do i=1,NX+1
       do j=1,NY+1
          if(j==1)then
             EY(i,j)=-(-3*phi(i,j)+4*phi(i,j+1)-phi(i,j+2))/(2*dy)
          else if(j==NY+1)then
             EY(i,j)=-(3*phi(i,j)-4*phi(i,j-1)+phi(i,j-2))/(2*dy)
          else
             EY(i,j)=-(phi(i,j+1)-phi(i,j-1))/(2*dx)
          end if
       end do
    end do



    !脳死差分法
    !do j=1,NY+1
    !   do i=1,NX+1
    !      EX(i,j)=-(phi(i+1,j)-phi(i,j))/dx
    !   end do
    !end do
    !do i=1,NX+1
    !   do j=1,NY+1
    !      EY(i,j)=-(phi(i,j+1)-phi(i,j))/dy
    !   end do
    !end do




    !!重み(Si/S)を用いて格子→粒子位置へ電場を与える
    !初期分布粒子に対して付与する
    do i=1,NX
       do j=1,NY
          do l=1,N_super_number
             hanteix=location_x(l)!!粒子のx座標
             hanteiy=location_y(l)!!粒子のy座標
             
             if(dx*(i-1)<=hanteix.and.hanteix<=dx*(i-1)+dx&
                  .and.dy*(j-1)<=hanteiy.and.hanteiy<=dy*(j-1)+dy)then
                S4=(hanteix-dx*(i-1))*(hanteiy-dy*(j-1))
                S3=(dx*i-hanteix)*(hanteiy-dy*(j-1))
                S2=(hanteix-dx*(i-1))*(dy*j-hanteiy)
                S1=(dx*i-hanteix)*(dy*j-hanteiy)
                !!対角の面積比から、電場を分配。
                EX_particle(l)=EX_particle(l)+(S1*EX(i,j)+S2*EX(i+1,j)+&
                     S3*EX(i,j+1)+S4*EX(i+1,j+1))/dx/dy
                EY_particle(l)=EY_particle(l)+(S1*EY(i,j)+S2*EY(i+1,j)+&
                     S3*EY(i,j+1)+S4*EY(i+1,j+1))/dx/dy
             end if
          end do
       end do
    end do


    !流入粒子に対して付与する
    do i=1,NX
       do j=1,NY
          do l=1,4*delta_Nall
             do t=1,tstep_count              
                hanteix=flow_loc_x(l,t)!!粒子のx座標
                hanteiy=flow_loc_y(l,t)!!粒子のy座標             
                if(dx*(i-1)<=hanteix.and.hanteix<=dx*(i-1)+dx&
                     .and.dy*(j-1)<=hanteiy.and.hanteiy<=dy*(j-1)+dy)then
                   S4=(hanteix-dx*(i-1))*(hanteiy-dy*(j-1))
                   S3=(dx*i-hanteix)*(hanteiy-dy*(j-1))
                   S2=(hanteix-dx*(i-1))*(dy*j-hanteiy)
                   S1=(dx*i-hanteix)*(dy*j-hanteiy)
                   !!対角の面積比から、電場を分配。
                   flow_EXP(l,t)=flow_EXP(l,t)+(S1*EX(i,j)+S2*EX(i+1,j)+&
                        S3*EX(i,j+1)+S4*EX(i+1,j+1))/dx/dy                 
                   flow_EYP(l,t)=flow_EYP(l,t)+(S1*EY(i,j)+S2*EY(i+1,j)+&
                        S3*EY(i,j+1)+S4*EY(i+1,j+1))/dx/dy
                end if
             end do
          end do
       end do
    end do


    
  end subroutine E_field



  subroutine E_field_output
    real(8)::max_Efield,currentE
    character(len=30)::fname_vectorE
    character(len=50)::fname_vectorE_part

    max_Efield=0.0d0
    do i=1,NX+1
       do j=1,NY+1
          currentE=sqrt(EX(i,j)**2+EY(i,j)**2)
          if(max_Efield<currentE)then
             max_Efield=currentE
          end if        
       end do
    end do
    

    write(fname_vectorE,'(a,I7,a)') &
         'vectorE_t=',nint(3*dble(tstep_count)),'[ns].dat' 
    open(100,file=fname_vectorE , status='replace')
    do i=1,NX+1
       do j=1,NY+1
          write(100,'(2F10.3,3E20.10)') dx*(i-1)*x_normal*1.0d2&
               ,dy*(j-1)*x_normal*1.0d2&
               ,EX(i,j)*2.0d0/max_Efield&
               ,EY(i,j)*2.0d0/max_Efield&
               ,sqrt((EX(i,j)*E_normal)**2+(EY(i,j)*E_normal)**2)
          !,EX(i,j)/max_Efield,EY(i,j)/max_Efield
          !,EX(i,j)*E_normal,EY(i,j)*E_normal
       end do
       write(100,*) 
    end do
    close(100)


    !write(fname_vectorE_part,'(a,I7,a)') &
    !     'vectorE_paticle_t=',nint(3*dble(tstep_count)),'[ns].dat'
    !open(200,file=fname_vectorE_part , status='replace')
    !do l=1,N_super_number
    !   if(lg_array_charge(l))then
    !      write(200,'(2E15.5)') EX_particle(l),EY_particle(l)
    !   end if
    !end do
    !close(200)


  end subroutine E_field_output
  !????????????????????????????????
  !????????????????????????????????
 
  
  subroutine Newton_eulerlow
    real(8)::ion_deltaloc_max,deltaloc_ion
    real(8)::ele_deltaloc_max,deltaloc_ele
    real(8)::ion_maxvel,ion_vel
    !!速度計算
    !!電子の場合(-Q0_super,SM0_super)
    !!イオンの場合(+Q0_super,SMXe_super)    do t=1,tstep_count
    !初期分布粒子
    deltaloc_ion=0.0d0
    ion_deltaloc_Max=0.0d0

    deltaloc_ele=0.0d0
    ele_deltaloc_max=0.0d0
    

    ion_vel=0.0d0
    ion_maxvel=0.0d0
    
    do l=1,N_super_number
       if(lg_array_charge(l))then
          velocity_x(l)=velocity_x(l)+(Q0_super/SMXe_super)*EX_particle(l)*dt
          velocity_y(l)=velocity_y(l)+(Q0_super/SMXe_super)*EY_particle(l)*dt

          ion_vel=sqrt((velocity_x(l))**2+(velocity_y(l))**2)
          if(ion_maxvel<ion_vel)then
             ion_maxvel=ion_vel
          end if
          
       end if
       
       if(.not.lg_array_charge(l))then         
          velocity_x(l)=velocity_x(l)+(-Q0_super/SM0_super)*EX_particle(l)*dt
          velocity_y(l)=velocity_y(l)+(-Q0_super/SM0_super)*EY_particle(l)*dt
          
       end if
    end do

    
    !!位置計算
    do l=1,N_super_number
       if(lg_array_charge(l))then
          location_x(l)=location_x(l)+velocity_x(l)*dt
          location_y(l)=location_y(l)+velocity_y(l)*dt

          deltaloc_ion=sqrt((velocity_x(l)*dt)**2+(velocity_y(l)*dt)**2)
          if(ion_deltaloc_max<deltaloc_ion)then
             ion_deltaloc_max=deltaloc_ion
          end if
          
       else if(.not.lg_array_charge(l))then
          location_x(l)=location_x(l)+velocity_x(l)*dt
          location_y(l)=location_y(l)+velocity_y(l)*dt

          deltaloc_ele=sqrt((velocity_x(l)*dt)**2+(velocity_y(l)*dt)**2)
          if(ele_deltaloc_max<deltaloc_ele)then
             ele_deltaloc_max=deltaloc_ele
          end if
       end if 
    end do

    !ionの最大速度,ぼーむ速度
    print'("#ion max velocity [m/s]")'
    print'(E20.10)',ion_maxvel*v_normal
    !ionの最大変位
    print'("#ion max delta location [m]")'
    print'(E20.10)',ion_deltaloc_max*x_normal

    !電子の最大変位
    print'("#electron max delta location [m]")'
    print'(E20.10)',ele_deltaloc_max*x_normal
    
!############################################
!############################################
!############################################
    !流入粒子
    !速度
    do l=1,4*delta_Nall
       do t=1,tstep_count
          flow_vel_x(l,t)=flow_vel_x(l,t)+(-Q0_super/SM0_super)*flow_EXP(l,t)*dt
          flow_vel_y(l,t)=flow_vel_y(l,t)+(-Q0_super/SM0_super)*flow_EYP(l,t)*dt
       end do
    end do

    !位置
    do l=1,4*delta_Nall
       do t=1,tstep_count
          flow_loc_x(l,t)=flow_loc_x(l,t)+flow_vel_x(l,t)*dt
          flow_loc_y(l,t)=flow_loc_y(l,t)+flow_vel_y(l,t)*dt
       end do
    end do
    

  end subroutine Newton_Eulerlow




    
  subroutine location_output
    character(len=30)::true_fname,false_fname,fname

    !ion location
    write(true_fname,'(a,I7.0,a)') 't_true='&
         ,nint(3*dble(tstep_count)),'[ns].dat'
    open(10,file=true_fname , status='replace')
    
    do l=1,N_super_number
       if(lg_array_charge(l))then
          write(10,'(2F20.8)') location_x(l)*x_normal*1.0d2,&
               location_y(l)*x_normal*1.0d2
       end if
    end do
    close(10)

    !#########################################################
    !#########################################################
    !electron location
    write(false_fname,'(a,I7.0,a)') 't_false=',nint(3*dble(tstep_count)),'[ns].dat'
    open(20,file=false_fname , status='replace')
    
    do l=1,N_super_number
       if(.not.lg_array_charge(l))then
          write(20,'(2F20.8)') location_x(l)*x_normal*1.0d2,&
               location_y(l)*x_normal*1.0d2
       end if
    end do

    write(20,*)
    
    do l=1,4*delta_Nall
       do t=1,tstep_count
          write(20,'(2F20.8)') flow_loc_x(l,t)*x_normal*1.0d2,&
               flow_loc_y(l,t)*x_normal*1.0d2
       end do
    end do
    close(20)
    !#########################################################
    !#########################################################

    
  end subroutine location_output



  
  subroutine initializing_euler_amounts
    real(8)::hanteix,hanteiy
    
    !乱数生成るーちん
    integer::clock_time,seed_size
    integer,allocatable::seed(:)
    real(8),allocatable::RN(:)
    
    allocate(RN(N_super_number))
    
    call system_clock(clock_time)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed=clock_time
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(RN(:))

    EX_particle(:)=0.0d0
    EY_particle(:)=0.0d0
    
    !!計算領域から飛び出した粒子を初期状態へ移行
    do l=1,N_super_number
       
       !!変更する粒子の情報だけ書き換えればよい
       !lg_array_loc=.true. =粒子は領域内に存在
       if(0.0d0<=location_x(l).and.location_x(l)<=XX.and.&
            0.0d0<=location_y(l).and.location_y(l)<=YY)then
          lg_array_loc(l)=.true.
       end if
       
       !計算対象外へ
       if((i_begin-1)*dx<=location_x(l).and.location_x(l)<=(i_end-1)*dx.and.&
            (j_begin-1)*dy<=location_y(l).and.location_y(l)<=(j_end-1)*dy)then
          location_x(l)=-100.0d0
          location_y(l)=-100.0d0          
          velocity_x(l)=0
          velocity_y(l)=0
          lg_array_loc(l)=.false.
       end if

       if(location_x(l)<=0.0d0.and.XX<=location_x(l).and.&
            location_y(l)<=0.0d0.and.YY<=location_y(l))then
          location_x(l)=-100.0d0
          location_y(l)=-100.0d0          
          velocity_x(l)=0
          velocity_y(l)=0
          lg_array_loc(l)=.false.
       end if
    end do


    flow_EXP(l,t)=0.0d0
    flow_EYP(l,t)=0.0d0
    
    do l=1,delta_Nall
       do t=1,tstep_count
          !!変更する粒子の情報だけ書き換えればよい
          !!!lg_array_loc=.true. =粒子は領域内に存在
          if(0.0d0<=flow_loc_x(l,t).and.flow_loc_x(l,t)<=XX.and.&
               0.0d0<=flow_loc_y(l,t).and.flow_loc_y(l,t)<=YY)then
             flow_lg_loc(l,t)=.true.
          end if
          !計算対象外へ
          if((i_begin-1)*dx<=flow_loc_x(l,t).and.flow_loc_x(l,t)<=(i_end-1)*dx.and.&
               (j_begin-1)*dy<=flow_loc_y(l,t).and.flow_loc_y(l,t)<=(j_end-1)*dy)then
             flow_loc_x(l,t)=-100.0d0
             flow_loc_y(l,t)=-100.0d0  
             flow_vel_x(l,t)=0.0d0
             flow_vel_y(l,t)=0.0d0
             flow_lg_loc(l,t)=.false.
          end if
          
          if(flow_loc_x(l,t)<=0.0d0.and.XX<=flow_loc_x(l,t).and.&
            flow_loc_y(l,t)<=0.0d0.and.YY<=flow_loc_y(l,t))then
             flow_loc_x(l,t)=-100.0d0
             flow_loc_y(l,t)=-100.0d0  
             flow_vel_x(l,t)=0.0d0
             flow_vel_y(l,t)=0.0d0
             flow_lg_loc(l,t)=.false.
          end if
       end do
    end do

       
       !薄膜へ吸収後は、速度の大きい成分の境界端から流入
       !if((i_begin-1)*dx<=location_x(l).and.location_x(l)<=(i_end-1)*dx.and.&
       !     (j_begin-1)*dy<=location_y(l).and.location_y(l)<=(j_end-1)*dy)then

          !vxがvyよりも大きければ、x軸から流入
          !if(abs(velocity_x(l))>abs(velocity_y(l)))then
             !!vxの正負で流入方向を決定
          !   if(velocity_x(l)>0.0d0)then
          !      location_x(l)=0.0d0
          !      location_y(l)=RN(l)*YY
          !   else
          !      location_x(l)=XX
          !      location_y(l)=RN(l)*YY
          !   end if
             
          !vyのほうが大きければy軸から流入
          !else
             !!vyの正負で流入方向を決定
          !   if(velocity_y(l)>0.0d0)then
          !      location_x(l)=RN(l)*XX
          !      location_y(l)=0.0d0
          !   else
          !      location_x(l)=RN(l)*XX
          !      location_y(l)=YY
          !   end if
          !end if
    
             
!周期境界条件          
!####################
       !領域から飛び出た粒子の速度を維持したまま
       !逆の境界端から流入
       !else if(location_x(l)<0.0d0)then
       !   location_x(l)=location_x(l)+XX
       !else if(XX<location_x(l))then
       !   location_x(l)=location_x(l)-XX
       !else if(location_y(l)<0.0d0)then
       !   location_y(l)=location_y(l)+YY
       !else if(YY<location_y(l))then
       !   location_y(l)=location_y(l)-YY
       !end if
!#####################       
       
          

       !境界から流出条件
       !if(location_x(l)<=0.0d0.or.XX<=location_x(l).or.&
       !     location_y(l)<=0.0d0.or.YY<=location_y(l))then
          !lg_array_loc=.false. =計算領域から飛び出た粒子
          !計算対象外へ
          !location_x(l)=-100.0d0
          !location_y(l)=-100.0d0          
          !velocity_x(l)=0
          !velocity_y(l)=0
       !   lg_array_loc(l)=.false.
       !end if


!############################
       !計算領域端反射条件 
       !if(location_x(l)<=0.0d0)then
       !   lg_array_loc(l)=.true.
          
       !   location_x(l)=-location_x(l)
       !   velocity_x(l)=-velocity_x(l)
       !end if
       
       !if(XX<=location_x(l))then
       !   lg_array_loc(l)=.true.

       !   location_x(l)=2*XX-location_x(l)
       !   velocity_x(l)=-velocity_x(l)
       !end if
       
       !if(location_y(l)<=0.0d0)then
       !   lg_array_loc(l)=.true.

       !   location_y(l)=-location_y(l)
       !   velocity_y(l)=-velocity_y(l)
       !end if
       
       !if(YY<=location_y(l))then
       !   lg_array_loc(l)=.true.

       !   location_y(l)=2*YY-location_y(l)
       !   velocity_y(l)=-velocity_y(l)
       !end if
!############################

    

          
          
  end subroutine initializing_euler_amounts















































  






































  









































  
  



  
  
  subroutine maxwellian_ele
    
    !number=N_super_number/2が計算領域の全超粒子数にあたる
    integer::number

    integer::NV      !step数 
    real(8)::dv=1.0d3 
    
    !dv毎の速度
    real(8),allocatable::maxwell_v(:)

    
    !速度格納配列
    real(8),allocatable::vx(:),vy(:),vx_r(:),vy_r(:)

    integer::vv
    integer::counter
    real(8)::tmp

    !すけーりんぐ予想
    real(8)::fv_max_1D_stat,vmax
    !定数部分
    real(8)::fv_const
    
                      !確率密度関数   累積確率密度関数
    real(8),allocatable::fv_max_1D(:),FV_max_cumu_1D(:)
    
    !###
    !任意の物理量＊f(v)の積分値
    real(8),allocatable::integral_vfv_none_v(:)
    !積分値 [-∞~∞]範囲の積分値
    real(8)::integral_vfv_none_finite
    !任意の物理量関数の、累積分布関数F(v)
    !F(v)=integral_vfv_none_v(:)/integral_vfv_none_finite
    real(8),allocatable::expect_Fv(:)
    !###

    !dv速度帯における粒子数counter
    real(8),allocatable::maxwell_v_count(:)
    
    !乱数生成routine
    real(8),allocatable::RN(:)
    integer,allocatable::int_RN(:)
    real(8)::R
    integer::clock_time,seed_size
    integer,allocatable::seed(:)
    
    character(len=50)::max_1d_fname,made_1d_fname


    fv_const=sqrt(1/pi)/vth_p3d
    number=N_super_number/2

    print'("Generating maxellian now (super electron)")'

    !すけーりんぐ予想
    !maxwellianがどの範囲?<v<?まで存在するか
    vmax=0
    do 
       fv_max_1D_stat=fv_const*exp(-vmax**2/(vth_p3d**2))
       vmax=vmax+dv
       !確率許容値
       if(fv_max_1D_stat<1.0d-10)exit
    end do
    vmax=vmax-dv
    !速度分布なので折り返した値 -vmax が最小値

    !step回数取得
    NV=nint(vmax/dv)
    
    
    allocate(maxwell_v(2*NV+2),fv_max_1D(2*NV+2),FV_max_cumu_1D(2*NV+2),&
         integral_vfv_none_v(2*NV+2),expect_Fv(2*NV+2))
    allocate(maxwell_v_count(2*NV+2))
    
    !!maxwell分布におけるv,f(v),F(v)を求める
    maxwell_v(:)=0
    fv_max_1D(:)=0
    FV_max_cumu_1D(:)=0
    integral_vfv_none_v(:)=0
    integral_vfv_none_finite=0
    maxwell_v_count(:)=0
  
    do vv=2,2*NV+1
       maxwell_v(vv)=-vmax+dv*(vv-1)
       fv_max_1D(vv)=fv_const*exp(-maxwell_v(vv)**2/(vth_p3d**2))
       FV_max_cumu_1D(vv)=FV_max_cumu_1D(vv-1)+fv_max_1D(vv)*dv
       !print'(2E15.5)',maxwell_v(vv),fv_max_1D(vv)*maxwell_v(vv)
    end do

    
    
    do vv=1,2*NV
       !積分値から、速度の累積分布関数F(v)を求める
       integral_vfv_none_v(vv)=integral_vfv_none_v(vv-1)&
            +abs(maxwell_v(vv)*fv_max_1D(vv)*dv)
    end do
    integral_vfv_none_finite=abs(integral_vfv_none_v(NV+1))/2

    !違う??????
    !delta_Nsuper=int(dt*dx*dy*integral_vfv_none_finite&
    !     *(ne2d**(1.5d0))/N_super/2)


    
    delta_Nsuper=int(dt*dx*integral_vfv_none_finite&
         *ne2d/N_super)
    delta_Nsuper_float=dble(dt*dx*integral_vfv_none_finite&
         *ne2d/N_super)

    print'("integral_vfv_-∞to∞",E15.5)', integral_vfv_none_finite
    print'("実際に使用するのは、 delta_Nsuper_float")'
    print'("delta_Nsuper(int) , delta_Nsuper_float(real)")'        
    print'(I15.5,E15.5)',delta_Nsuper,delta_Nsuper_float
    delta_Nall=nint(delta_Nsuper_float*NX)
    print'("1つの境界から流入させる総電子数")'
    print'("delta_Nall=nint(delta_Nsuper_float*NX)")'
    print'("delta_Nall",I15.5)',delta_Nall
    
    
    !F(v)計算
    expect_Fv(:)=0
    expect_Fv(:)=integral_vfv_none_v(:)/integral_vfv_none_finite
    !do vv=1,2*NV
    !   expect_Fv(vv)=integral_vfv_none_v(vv)/integral_vfv_none_finite
    !   !print'(2E15.5)',maxwell_v(vv),expect_Fv(vv)
    !end do
    









    
    !###################################
    !maxwelian生成routine
    allocate(RN(number))
    call system_clock(clock_time)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed=clock_time
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(RN(:))


    maxwell_v_count(:)=0
    !速度帯dvの間の存在数をかうんと
    do l=1,number
       do vv=1,2*NV
          !if(expect_Fv(vv)<RN(l).and.RN(l)<expect_Fv(vv+1))then
          if(FV_max_cumu_1D(vv)<RN(l).and.RN(l)<FV_max_cumu_1D(vv+1))then
             maxwell_v_count(vv)=maxwell_v_count(vv)+1
          end if
       end do
    end do


!###################################################
    !1Dのmaxwellian
    !理論値 theoretical values
    write(max_1d_fname,'(a)') 'ini_fv_false_1d_theory.dat'
    open(10,file=max_1d_fname , status='replace')
    do vv=2,2*NV+1
       write(10,'(2E20.8)') maxwell_v(vv), &
            fv_max_1D(vv)!abs(maxwell_v(vv)*fv_max_1D(vv))/vth_p3d!/integral_vfv_none_finite) 
    end do
    close(10)
    
    
    !do i=1,50
    i=1
    !作成したmaxwell分布関数 1D
       write(made_1d_fname,'(a,I7.0,a)') 'ini_fv_false_1d_made',i,'.dat'
       open(100,file=made_1d_fname , status='replace')
       do vv=2,2*NV+1
          write(100,'(2E20.8)') maxwell_v(vv), &
               dble(maxwell_v_count(vv))/dble(number)/dv!/vth_p3d/dv*3.2
       end do
       close(100)
    !end do
    
!#########################################

    
    !countした粒子ぶんだけ、その速度帯に代入
    allocate(vx(number),vx_r(number))
    counter=0
    do vv=1,2*NV
       if(.not.maxwell_v_count(vv)==0)then
          do while(.not.maxwell_v_count(vv)==0)
             counter=counter+1
             vx(counter)=maxwell_v(vv)+0.5d0*dv
             maxwell_v_count(vv)=maxwell_v_count(vv)-1
          end do
       end if
    end do
    counter=1
    maxwell_v_count(:)=0


    call system_clock(clock_time)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed=clock_time
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(RN(:))
    
    !Fisher-Yates Shuffle
    !しゃっふる あるごりずむ
    !しゃっふる後の格納配列を用意する事で、poiner error回避
    do i=1,number
       j=mod(RN(i)*dble(number),dble(number))      
       vx_r(i)=vx(j)   
    end do
    !########ここまで


    !同じるーちんで今度はvyを求める
    call system_clock(clock_time)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed=clock_time
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(RN(:))


    maxwell_v_count(:)=0
    !速度帯dvの間の存在数をかうんと
    do l=1,number
       do vv=1,2*NV
          !if(expect_Fv(vv)<RN(l).and.RN(l)<expect_Fv(vv+1))then
          if(FV_max_cumu_1D(vv)<RN(l).and.RN(l)<FV_max_cumu_1D(vv+1))then
             maxwell_v_count(vv)=maxwell_v_count(vv)+1
          end if
       end do
    end do


    !countした粒子ぶんだけ、その速度帯に代入
    allocate(vy(number),vy_r(number))
    counter=0
    do vv=1,2*NV
       if(.not.maxwell_v_count(vv)==0)then
          do while(.not.maxwell_v_count(vv)==0)
             counter=counter+1
             vy(counter)=maxwell_v(vv)+0.5d0*dv
             maxwell_v_count(vv)=maxwell_v_count(vv)-1
          end do
       end if
    end do
    counter=1
    maxwell_v_count(:)=0


    call system_clock(clock_time)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed=clock_time
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(RN(:))

    
    !Fisher-Yates Shuffle
    !しゃっふる あるごりずむ
    !しゃっふる後の格納配列を用意する事で、poiner error回避
    do i=1,number
       j=mod(RN(i)*dble(number),dble(number))      
       vy_r(i)=vy(j)   
    end do

    
    !共用変数へ代入
    counter=1
    do l=1,N_super_number
       !負電荷判定
       if(.not.lg_array_charge(l))then
          velocity_x(l)=vx_r(counter)
          velocity_y(l)=vy_r(counter)
          counter=counter+1
       end if
    end do










    

!###################################################
    !ionにも 代入してみる
    !counter=1
    !do l=1,N_super_number
       !正電荷判定
    !   if(lg_array_charge(l))then
    !      velocity_x(l)=vx_r(counter)*1.0d-03
    !      velocity_y(l)=vy_r(counter)*1.0d-03
    !      counter=counter+1
    !   end if
    !end do
!###################################################

    !ionはcold
    do l=1,N_super_number
       if(lg_array_charge(l))then
          velocity_x(l)=0.0d0
          velocity_y(l)=0.0d0
       end if
    end do

























    







    

    


    

    
    !###################
    !指定した速度分布関数から、流入速度を決定する
    call system_clock(clock_time)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed=clock_time
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(RN(:))
    allocate(vx_inflow_neg(delta_Nall),vx_inflow_pos(delta_Nall),&
         vy_inflow_neg(delta_Nall),vy_inflow_pos(delta_Nall))
    allocate(v_inflow_thermal(delta_Nall))

    
    

    !指定した速度分布関数から流入速度を求める
    !###############################################    

    !                        image 
    !                  (vth,vy_inflow_pos)
    !                    ###############
    !                    #             #
    !                    #             #
    !(vx_inflow_neg,vth) #             #(vx_inflow_pos,vth)
    !                    #             #
    !                    #             #
    !                    #             #
    !                    ###############
    !                  (vth,vy_inflow_neg)
    !
    
   
    do i=1,delta_Nall
       j=mod(RN(i)*dble(number/2),dble(number/2))
       vx_inflow_pos(i)=vx(j)
    end do

    do i=1,delta_Nall
       j=mod(RN(i)*dble(number/2),dble(number/2))
       vx_inflow_neg(i)=vx(j+(number/2))
    end do
    
    
    call system_clock(clock_time)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed=clock_time
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(RN(:))
    
    !同様にvy
    do i=1,delta_Nall
       j=mod(RN(i)*dble(number/2),dble(number/2))
       vy_inflow_pos(i)=vy(j)
    end do

    do i=1,delta_Nall
       j=mod(RN(i)*dble(number/2),dble(number/2))
       vy_inflow_neg(i)=vy(j+(number/2))
    end do

    !同様にvth
    do i=1,delta_Nall
       j=mod(RN(i)*dble(number),dble(number))
       v_inflow_thermal(i)=vy(j)
          !どちらも同じ分布関数なので、vxでもvyでも構わない
    end do
   
    

    !print out
    !print'("conb:(negx,y)")'
    !do l=1,delta_Nall
    !   print'(2E15.5)',vx_inflow_neg(l),v_inflow_thermal(l)
    !end do
    
    !print'("conb:(posx,y)")'
    !do l=1,delta_Nall
    !   print'(2E15.5)',vx_inflow_pos(l),v_inflow_thermal(l)
    !end do
    
    !print'("conb:(x,negy)")'
    !do l=1,delta_Nall
    !   print'(2E15.5)',v_inflow_thermal(l),vy_inflow_neg(l)
    !end do
    
    !print'("conb:(x,posy)")'
    !do l=1,delta_Nall
    !   print'(2E15.5)',v_inflow_thermal(l),vy_inflow_pos(l)
    !end do
!###############################################    

    
  end subroutine maxwellian_ele























  




    !fiting対象表示用
  subroutine maxwellian_fv2d
    character(len=50)::fv_false_first_fname
    real(8)::fv_max_2d_stat,fv_const2d,vmax
    integer::vv,NV
    real(8)::dv=1.0d3
    real(8),allocatable::maxwell_v(:),fv_max_2D(:),FV_max_cumu_2D(:)
    real(8),allocatable::integral_vfv_none_v(:),expect_Fv(:),maxwell_v_count(:)
    real(8)::integral_vfv_none_finite


    fv_const2d=2/(vth_p3d**2)
    !2d|v|のmaxwellianに代入した electron f(v)-vを求める
    !すけーりんぐ予想
    !maxwellianがどの範囲?<v<?まで存在するか

    vmax=dv
    do
       fv_max_2D_stat=fv_const2d*vmax*exp(-vmax**2/vth_p3d**2)
       vmax=vmax+dv
       !確率許容値
       if(fv_max_2D_stat<1.0d-10)exit
    end do
    vmax=vmax-dv

    !step回数取得
    NV=nint(vmax/dv)
    
    
    allocate(maxwell_v(NV+2),fv_max_2D(NV+2),FV_max_cumu_2D(NV+2),&
         integral_vfv_none_v(NV+2),expect_Fv(NV+2))
    allocate(maxwell_v_count(NV+2))
    
    !!maxwell分布におけるv,f(v),F(v)を求める
    maxwell_v(:)=0
    fv_max_2D(:)=0
    FV_max_cumu_2D(:)=0
    integral_vfv_none_v(:)=0
    integral_vfv_none_finite=0
    maxwell_v_count(:)=0


    do vv=1,NV+1
       maxwell_v(vv)=dv*(vv-1)
       fv_max_2D(vv)=fv_const2d*maxwell_v(vv)*exp(-maxwell_v(vv)**2/vth_p3d**2)
       FV_max_cumu_2D(vv)=FV_max_cumu_2D(vv-1)+fv_max_2D(vv)*dv
    end do

    !積分値
    do vv=2,NV+1
       !積分値から、期待値の累積分布関数F(v)を求める
       integral_vfv_none_v(vv)=integral_vfv_none_v(vv-1)&
            +abs(maxwell_v(vv)*fv_max_2D(vv)*dv)
    end do
    integral_vfv_none_finite=abs(integral_vfv_none_v(NV+1)) 


    !速度の2D確率分布関数 理論値
    !ini_fv_false_2d_theory.dat
    write(fv_false_first_fname,'(a)') 'ini_fv_false_2d_theory.dat'
    open(100,file=fv_false_first_fname , status='replace')
    do vv=1,NV+1
       write(100,'(2E20.8)') maxwell_v(vv),&
            fv_max_2D(vv)!*maxwell_v(vv)/integral_vfv_none_finite
    end do
    close(100)

  end subroutine maxwellian_fv2d
  














  subroutine ini_ele_fv_output
    !初期分布表示用
    character(len=50)::fv_false_fname,vx_buffer_fname,vy_buffer_fname
    integer::vv,NV,number,counter
    integer::tstep=0
    real(8)::dv=1.0d3
    real(8)::vmax
    real(8)::fv_max_2D,fv_const2D
    real(8),allocatable::v(:)
    integer,allocatable::maxwell_v_count(:)
    real(8),allocatable::integral_vfv(:)
    real(8)::integral_vfv_finite

    !乱数生成るーちん
    integer::clock_time,seed_size
    integer,allocatable::seed(:)
    real(8),allocatable::RN(:)
    fv_const2d=2/(vth_p3d**2)


    !すけーりんぐ予想
    !maxwellianがどの範囲?<v<?まで存在するか
    vmax=dv
    do
       fv_max_2D=fv_const2d*vmax*exp(-vmax**2/vth_p3d**2)
       vmax=vmax+dv
       !確率許容値
       if(fv_max_2D<1.0d-10)exit
    end do
    vmax=vmax-dv

    NV=nint(vmax/dv)
    allocate(maxwell_v_count(NV+2),integral_vfv(NV+2))
    maxwell_v_count(:)=0
    integral_vfv(:)=0

    !境界条件含めた粒子数をcount
    counter=0
    do l=1,N_super_number
       if(.not.lg_array_charge(l))then
          counter=counter+1
       end if    
    end do
    
    number=counter
    allocate(v(number+2))
    v(:)=0
 

    allocate(RN(N_super_number))
    call system_clock(clock_time)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    call random_seed(put=seed+tstep)
    deallocate(seed)
    call random_number(RN(:))


    
    counter=1
    do l=1,N_super_number
       if(.not.lg_array_charge(l).and.RN(l)>0.1d0)then
          !粒子の速さ格納配列
          v(counter)=sqrt((velocity_x(l))**2+&
               (velocity_y(l))**2)
          counter=counter+1
       end if    
    end do

    
    do l=1,number
       do vv=1,NV         
          if(dv*(vv-1)<=v(l).and.v(l)<=dv*vv)then
             maxwell_v_count(vv)=maxwell_v_count(vv)+1
          end if
       end do
    end do


    do vv=2,NV
       !積分値から、期待値の累積分布関数F(v)を求める
       integral_vfv(vv)=integral_vfv(vv-1)+&
            dble(abs(dv*(vv-1)*dble(maxwell_v_count(vv))/dble(number)))    
    end do
    integral_vfv_finite=integral_vfv(NV)



    tstep=tstep+1
    !file output
    !作成した2D maxwell分布 @t=0
    !任意の回数出力(回数はmain関数を参照)
    !ini_fv_false_2d_made_数字.dat 
    write(fv_false_fname,'(a,I7.0,a)')'ini_fv_false_2d_made' ,tstep,'.dat'
    open(100,file=fv_false_fname , status='replace')
    do vv=2,NV
       write(100,'(2E20.8)') dv*(vv-1), &
            dble(maxwell_v_count(vv))/dble(number)/dv!/integral_vfv_finite!
    end do
    close(100)

    !file output
    !作成した速度をそのまま出力
    !vxbuffer  数字.dat  
    write(vx_buffer_fname,'(a,I7.0,a)') 'vxbuffer' ,tstep,'.dat'
    open(10,file=vx_buffer_fname , status='replace')
    do l=1,N_super_number
       if(.not.lg_array_charge(l).and.RN(l)>0.0d0)then
          write(10,'(2F20.8)') velocity_x(l)
       end if
    end do
    close(10)

    !file output
    !vybuffer  数字.dat
    !正味いらない
    !write(vy_buffer_fname,'(a,I7.0,a)') 'vybuffer' ,tstep,'.dat'
    !open(20,file=vy_buffer_fname , status='replace')
    !do l=1,N_super_number
    !   if(.not.lg_array_charge(l).and.RN(l)>0.7d0)then
    !      write(10,'(2F20.8)') velocity_y(l)
    !   end if
    !end do
    !close(20)
  end subroutine ini_ele_fv_output



  
  






  subroutine C_density_distribute_output
    character(len=30)::fname_Cden

    
    write(fname_Cden,'(a,I7,a)') 'Cden=',nint(3*dble(tstep_count)),'[ns].csv'
    
    open(100,file=fname_Cden , status='replace')

    do i=1,NX+1
       do j=1,NY+1
          write(100,'(E15.3,",",E15.3,",",E20.10)')&
               dx*(i-1)*x_normal*1.0d2,&
               dy*(j-1)*x_normal*1.0d2,&
               charge(i,j)*rho_normal
       end do
       write(100,*) !xの値が変更する際、改行しないとgnuplotは正しく読み込んでくれない。
    end do
    close(100)
  end subroutine C_density_distribute_output
  





  







  
end module subprog



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
  use subprog
  implicit none
  !mainでの総時間めっしゅ,るーぷ用変数
  integer::NT_main,tmain
  !file output span
  integer::file_out_span_main
  !初期速度分布抽出回数
  integer::fvout
  !mainでの電子温度,電子質量,Xeいおん質量
  real(8)::Te_main,SM0_super_main,SMXe_super_main


  
  call input_allocate !!入力&動的割付
  call passing_to_main(NT_main,file_out_span_main&
       ,Te_main,SM0_super_main,SMXe_super_main)



  !電荷の正負決定
  call charge_decided
  
  !maxwellian生成(超粒子電子)
  call maxwellian_ele

!######################################
  !fiting対象関数出力
  !call maxwellian_fv2d

  
  !初期分布から粒子をらんだむに抽出後、速度分布作成
  !任意回数出力
  !do fvout=1,50
  !   call ini_ele_fv_output
  !end do
!###########################################3


  !物理量の規格化
  call normalization

  !初期一様分布粒子
  !初期条件作成
  call initial_conditions
     
  do tmain=1,NT_main
     call counter  !counter管理
     call initial_conditions_inflow !!流入粒子の初期条件
     call C_density_distribute !!電荷密度を、重みを用いてめっしゅに与える。
     call poisson_Gauss_Seidel !!poisson方程式の差分化 (GaussSeidel法


     call E_field !!重みを用いてめっしゅ→粒子位置へ与える
     
     call Newton_eulerlow !!euler法により、運動方程式を解き、粒子の位置,速度を求める
     

     !##dt=3.0[ns]*数値stepで出力
     if(mod(tmain,file_out_span_main)==0)then
        call C_density_distribute_output
        call E_field_output
        call phidis_output !!電位分布表示 
        call location_output !!粒子位置表示
        
        print'()'
        print'()'
     end if
      
     
   
     call initializing_euler_amounts !変数の初期化


     !時間経過の速度分布出力
     !call ele_fv_output

     
     print'()'
  end do
  print'("計算を終了します")'
end program main

