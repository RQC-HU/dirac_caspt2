! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE print_head_casci

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables, only: rank

    Implicit NONE

    if (rank == 0) then
        print *, ''
        print '(84A)', '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        print '(84A)', '@        ____  _____ _        _  _____ _____     _____ ____ _____ ___ ____         @'
        print '(84A)', '@       |  _ \| ____| |      / \|_   _|_ _\ \   / /_ _/ ___|_   _|_ _/ ___|        @'
        print '(84A)', '@       | |_) |  _| | |     / _ \ | |  | | \ \ / / | |\___ \ | |  | | |            @'
        print '(84A)', '@       |  _ <| |___| |___ / ___ \| |  | |  \ V /  | | ___) || |  | | |___         @'
        print '(84A)', '@       |_|_\_\_____|_____/_/___\_\_|_|___|_ \_/  |___|____/ |_| |___\____|        @'
        print '(84A)', '@                         ____    _    ____   ____ ___                             @'
        print '(84A)', '@                        / ___|  / \  / ___| / ___|_ _|                            @'
        print '(84A)', '@                       | |     / _ \ \___ \| |    | |                             @'
        print '(84A)', '@                       | |___ / ___ \ ___) | |___ | |                             @'
        print '(84A)', '@                        \____/_/   \_\____/ \____|___|                            @'
        print '(84A)', '@                                                                                  @'
        print '(84A)', '@       Developped by Kohei Noda, Yasuto Masuda, Sumika Iwamuro, Minori Abe        @'
        print '(84A)', '@               Hiroshima University & Tokyo Metropolitan Univeristy               @'
        print '(84A)', '@                     https://github.com/RQC-HU/dirac_caspt2                       @'
        print '(84A)', '@                                                                                  @'
        print '(84A)', '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        print '(84A)', '@                                                                                  @'
        print '(84A)', '@                MMMN,                                    .MMp                     @'
        print '(84A)', '@               .#  ?YMa,                                dM\Tb                     @'
        print '(84A)', '@               .F    ,MMm,                             dM#  .W,                   @'
        print '(84A)', '@               (F      TMMp                           JMMt    7,                  @'
        print '(84A)', '@               J]       4MMp                      .dMM@       M[                  @'
        print '(84A)', '@               M\        ?MMh.                    .(MMMF       -N.                @'
        print '(84A)', '@              .M          qMMMN+.............(+gMMMMM@^ (((()  MN                 @'
        print '(84A)', '@              (F          .MT"^ _~?77?!~~???77"^ ``,!           JM!               @'
        print '(84A)', '@              M]                                                  TN,             @'
        print '(84A)', '@             .M}                                                  VMp             @'
        print '(84A)', '@             MM        .gMMMNJ.                 ..MMMMMa,          UM2            @'
        print '(84A)', '@            JP       .MMMMMMMMMN,              .M^MMMMN,TN,         UMp           @'
        print '(84A)', '@           Jb        M@.MMMMMM)?Mh           .M^ MMMMMM] ?h          MMb          @'
        print '(84A)', '@          .M^       .Hp 4MMMM!  .Mb          .M  (MMMM) ..9]      ...(MMNCASPT2   @'
        print '(84A)', '@         .Mh.(J.,     TMNMMMMNM"^ 7`         TWgJuMMNg&M"     .d"""777MM!~??`     @'
        print '(84A)', '@    ..gMMMM"^ .^          ?^                       ^^?                q@]         @'
        print '(84A)', '@  .MM#"` MF        ^^^               ......             ^^^           -b          @'
        print '(84A)', '@ , "     HN{                         WMMMM\   ..&,.           MNMNH96jM@          @'
        print '(84A)', '@         .Mh,                 .       /M@^        //              _7^WMMMMNg,.    @'
        print '(84A)', '@          dMugNMMMMMM9         .%      M:       ..M`                J@    74,     @'
        print '(84A)', '@        ..MMM&       ..        TMa....M9WNmQggN#^?          jQNgg+..MF            @'
        print '(84A)', '@      .M@^ .dMMo..gd ^`            `                          ?jM RASPT2.         @'
        print '(84A)', '@     dF      ~jMM#1.               .N,      .N=        `   ..g#      .7           @'
        print '(84A)', '@             .MMTTMMNm,.            -MMMWyM9-         ..+MMMT=)  dMM      mi      @'
        print '(84A)', '@           .dMF    7?W""9g...,                 ...J"=`y"D`         FFM            @'
        print '(84A)', '@          .MMF                 ^    """"""""   ^                      KK.         @'
        print '(84A)', '@                                                                                  @'
        print '(84A)', '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        print '(84A)', '@                                                                                  @'
        print '(84A)', '@       This program is intended to be used after HF and MO transformation         @'
        print '(84A)', '@       calculations uisng the DIRAC software.                                     @'
        print '(84A)', '@       https://www.diracprogram.org/doku.php?id=start                             @'
        print '(84A)', '@                                                                                  @'
        print '(84A)', '@       You can use, redistribute and/or modify this program under the terms       @'
        print '(84A)', '@       of the GNU Lesser General Public License version 2.1 as published          @'
        print '(84A)', '@       by the Free Software Foundation.                                           @'
        print '(84A)', '@                                                                                  @'
        print '(84A)', '@       Please site the following paper when you publish the data obtained         @'
        print '(84A)', '@       using this program.                                                        @'
        print '(84A)', '@       "xxx, xxx, xxx      "                                                      @'
        print '(84A)', '@                                                                                  @'
        print '(84A)', '@       We do not take any responsibility for the results produced                 @'
        print '(84A)', '@       by this program.                                                           @'
        print '(84A)', '@                                                                                  @'
        print '(84A)', '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        print *, ''
    end if

end subroutine print_head_casci
