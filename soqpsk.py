"""
soqpsk.py

Description: 

Donald MacIntyre - djm4912
"""

import math
import random

# This generates the bits of SUAS_DDL preamble
# in accordance with the SUAS DDL specification
def getSUAS_DDL_Preamble(bits):
    # RampUp - 1Byte
    # Guard - 4 Bytes
    # Sync1to30 - 4Bytes*30
    # Inverse Sync - 4Bytes

    BITS_PER_BYTE = 8
    # Rampup
    for i in range(BITS_PER_BYTE):
        bits.append(0)
    # Guard
    for i in range(4*BITS_PER_BYTE):
        bits.append(0)
    #SYNC1to30
    sync = [1,0,0,0,    # 8
            0,0,0,1,    # 1
            0,1,0,1,    # 5
            1,1,0,1,    # D
            1,0,0,0,    # 8
            0,1,1,1,    # 7
            1,0,0,1,    # 9
            0,0,1,1]    # 3

    for i in range(30):
        for j in sync:
            bits.append(j)
    # INVERSE SYNC
    inv_sync = []
    for i in sync:
        if i == 0:
            inv_sync.append(1)
        else:
            inv_sync.append(0)

    for i in inv_sync:
        bits.append(i)
    return bits

# This gets a SUAS payload from an
# input file. Note that the input file is expected to
# be in the format of 1 bit per line. With the data being
# RS encoded and scrambled
def getSUAS_Payload(bits):

    f = open('ddl_frame1_rs_encoded_and_scrambled.txt', 'r')

    for b in f:
        b = b.strip()
        bits.append(int(b))

    f.close()
    return bits

def getSUAS_DDL_Postamble(bits):
    # Guard - 4 Bytes
    # Ramp Down - 1 Byte
    BITS_PER_BYTE = 8
    #Guard
    for i in range(4*BITS_PER_BYTE):
        bits.append(0)

    #Ramp Down
    for i in range(BITS_PER_BYTE):
        bits.append(0)
    return bits

# Calculate the sinc function
# i.e. sin(pi*x)/pi*x
def sinc(x):
    if (math.pi*x) == 0:
        return 1
    else:
        return math.sin(math.pi*x)/(math.pi*x)

# This is the table used to turn the serial bit stream
# into a frequency impulse (alpha_k) for the shaping filter.
# This table is used during even bit intervals and is
# Table 30 from the SUAS-DDL specifiation
def evenBitInterval(i, qpast, ipast):
    if ((i  == 0) and (qpast == 0) and (ipast == 0)):
        return 0
    if ((i  == 0) and (qpast == 0) and (ipast == 1)):
        return -1
    if ((i  == 0) and (qpast == 1) and (ipast == 0)):
        return 0
    if ((i  == 0) and (qpast == 1) and (ipast == 1)):
        return 1
    if ((i  == 1) and (qpast == 0) and (ipast == 0)):
        return 1
    if ((i  == 1) and (qpast == 0) and (ipast == 1)):
        return 0
    if ((i  == 1) and (qpast == 1) and (ipast == 0)):
        return -1
    if ((i  == 1) and (qpast == 1) and (ipast == 1)):
        return 0
    
# This is the table used to turn the serial bit stream
# into a frequency impulse (alpha_k) for the shaping filter.
# This table is used during odd bit intervals and is
# Table 30 from the SUAS-DDL specifiation
def oddBitInterval(q, ipast, qpast):
    if ((q  == 0) and (ipast == 0) and (qpast == 0)):
        return 0
    if ((q  == 0) and (ipast == 0) and (qpast == 1)):
        return 1
    if ((q  == 0) and (ipast == 1) and (qpast == 0)):
        return 0
    if ((q  == 0) and (ipast == 1) and (qpast == 1)):
        return -1
    if ((q  == 1) and (ipast == 0) and (qpast == 0)):
        return -1
    if ((q  == 1) and (ipast == 0) and (qpast == 1)):
        return 0
    if ((q  == 1) and (ipast == 1) and (qpast == 0)):
        return 1
    if ((q  == 1) and (ipast == 1) and (qpast == 1)):
        return 0

# This generates the filters used in modulation. These filters define
# how the frequency impulses are shaped for transmit.
def generateModulationFilters(p, B, T1, T2, SAMPLES_PER_SYMBOL):

    # Generate the window for the shaping filter
    # This is the w(t) filter as specified in the SUAS DDL specification
    w = []

    # Windowing Filter is non-zero in value from -(T1+T2) to (T1+T2)
    # So in the case of SUAS runs over four SYMBOL intervals.
    # In creation of the windowing filter generate values for six symbols
    # instead of 4. This is not necessary and will result in values
    # of zero for the range outside of -(T1+T2) to (T1+T2)
    for i in range(-3*SAMPLES_PER_SYMBOL,3*SAMPLES_PER_SYMBOL+1):
        # This normalizes t... i.e. t in this script is now equal to
        # t/T as defined in the SUAS DDL specification
        t = i / SAMPLES_PER_SYMBOL
        # Region 1 - abs(t/T) < T1
        if abs(t) < T1:
            w.append(1)
        # Region 2 - T1 < abs(t/T) < T1 + T2
        elif abs(t) < (T1 + T2):
            w.append(0.5+0.5*math.cos(math.pi*(abs(t)-T1)/T2))
        # Region 3 - T1+T2 < abs(t/T)
        else:
            w.append(0)

    # Generate the pulse shaping filter. Note that this pulse shaping
    # filter will eventually be truncated by the windowed filter above.
    # At a high level frequency impulses will be filtered by this filter.
    # This convolution will generate a "phase-pulse".  This phase pulse
    # accumulated across a symbol is always -pi/2, 0 or pi/2. 
    n = []
    # A - This paramater controls the total phase shift that is achieved
    # from the shaping filter. This value must be set such that the total
    # phase shift for one frequency impulse is pi/2. Will account for this later
    # after the shaping filter has been windowed
    A = 1
    for i in range(-3*SAMPLES_PER_SYMBOL,3*SAMPLES_PER_SYMBOL+1):
        # This normalizes t... i.e. t in this script is now equal to
        # t/T as defined in the SUAS DDL specification
        t = i / SAMPLES_PER_SYMBOL
        # n(t) right from the SUAS DDL Specification
        n.append(((A*math.cos(math.pi*p*B*t))/(1-4*(p*B*t)**2))*(sinc(B*t)))

    # n(t) and w(t) have been generated. Now need to window n(t) with w(t)
    # to create g(t). g(t) is the final pulse shaping filter.
    # g(t) indicates at any given time what is the phase change at a particular
    # instance in time. 
    g = []
    # Do the windowing here of n(t) thus generating g(t)
    for i in range(len(n)):
        g.append(w[i]*n[i])
    # Earlier assumed A = 1
    # Need to noramilize so for one bit have a phase shift of pi/2 radians
    # i.e. the total phase change throughout all of g must be pi/2
    # Accumulate the total phase change in g
    total_phase = sum(g)
    # Determine a scale factor that can be used such that the total sum
    # of g(i.e. total_phase) can be normalized to sum to pi/2
    scale = (math.pi/2)/total_phase
    # Finally scale g by the aforementioned scaling factor.
    # g will now match the specification with (A) being set such that
    # accumulating over one frequency impulse generates a phase shift of
    # pi/2 radians
    for i in range(len(g)):
        g[i] *= scale

    # Generate the phase pulse.  The phase pulse is not used currently used
    # anywhere but is helpful in understanding what is going on.
    # The phase pulse is just an accumulation of phase change. So unlike
    # g which indicats the phase change at a specific time, the phase pulse
    # indicates the total phase at a given time
    phase_pulse = []
    last_phase = 0
    # Generate the phase pulse, Accumulation of all phase change seen
    # to this point + the phase change for this instance in time
    for i in range(len(g)):
        phase_pulse.append(g[i] + last_phase)
        last_phase = phase_pulse[i]

    # Returning the shaping filter and the phase_pulse
    return g, phase_pulse

# This function takes a serial bitstream and turns said bits
# into frequency impulses (alpha_k).
def getFrequencyImpulses(bits):
    alpha = []
    isEven = 1
    # Initalize lasti and lastq to zero for first time through the
    # modulation. This is valid as prior to sync, 00 is used for
    # rampup and guard sequence
    lasti = 0
    lastq = 0
    for j in range(len(bits)):
        # Grab the latest bit to process
        b = bits[j]
        # I bit Processing
        if isEven == 1:
            i = bits[j]
            # Determine the frequency impulse using even table
            alpha.append(evenBitInterval(i,lastq,lasti))
            # Move this i to last i
            lasti = i
        # Q bit Processing
        else:
            q = bits[j]
            # Determine the frequency impulse using odd table
            alpha.append(oddBitInterval(q,lasti,lastq))
            # Move this q to last q
            lastq = q
        # Toggle between even and odd bit intervals
        isEven ^= 1
    # Return the frequency impulses
    return alpha

# Normalization to keep angle in radians between -pi and pi
def normalizeAngleRadians(a):
    while (a < -math.pi):
        a += 2*math.pi
    while (a > math.pi):
        a -= 2*math.pi
    return a

# Function to perform convolution.
# No scaling and shifting... instead take impulse response and for each
# impulse in the input sequence (x), shift the impulse response to appropriate
# location and scale by input sequence. Convolution is summation of all impulse
# responses for each impulse in x
def conv(x,h):
    res = [0]*(len(x)+len(h)-1)
    for i in range(len(h)):
        for j in range(len(x)):
            res[i+j] += h[i]*x[j]
    return res       

def main() :

    # SUAS DDL Waveform Paramaters
    # Defined in SUAS DDL Specification
    p = 1
    B = 1.35
    T1 = 1.4
    T2 = 0.6

    # This is the number of samples per symbol (not samples per bit).
    # Note that a symbol consists of an I and Q bit.
    # Therefore the samples ber bit will be half of the samples per symbol.
    SAMPLES_PER_SYMBOL = 4

    # Generate phase shaping filter to be used in modulation
    g,phase_pulse = generateModulationFilters(p,B,T1,T2,SAMPLES_PER_SYMBOL)

    # Generate a random bit sequence of length NUM_BITS
    
    bits = []
    bits = getSUAS_DDL_Preamble(bits)

    bits = getSUAS_Payload(bits)
    
    # Add the SUAS posamble
    bits = getSUAS_DDL_Postamble(bits)
    
    # Turn the bits into frequency impulses
    alpha = getFrequencyImpulses(bits)

    # With the frequency impulses generated, next step is to convolve
    # aforementioned impulses with the phase shaping filter. Need to align
    # the frequency impulses such that the occur at a rate of the bit interval.
    # Therefore for each frequency impulse, need to zero pad with
    # (SAMPLES_PER_SYMBOL//2)-1 zeros. This will ensure that each
    # SYMBOL_INTERVAL will have two frequency pulses associated with it
    # one for I and one for Q. 
    alpha_impulses_one_per_symb = []
    for a in alpha:
        alpha_impulses_one_per_symb.append(a)
        for j in range((SAMPLES_PER_SYMBOL//2)-1):
            alpha_impulses_one_per_symb.append(0)
            
    # Generate the modulated phase. This is just the summation (convolution) of
    # all phase pulses for each individual frequency impulse.
    # Note that since the impulse response for an individual frequency impulse runs
    # over multiple symbol intervals and a new frequency impulse is generated
    # at a rate of the bit interval there is signficant inter-symbol-interfernce
    phase_delta = conv(alpha_impulses_one_per_symb,g)

    # phase_delta contains the phase change at a given point in time.
    # Would like to accumulate that phase change to be able to determine the
    # "phase" at any given point in time. At this point total_phase can be used
    # to generate the "traditional" SOQPSK constellation.
    total_phase = []
    # Adding in an offset of pi/4. SUAS DDL specification puts symbols at
    # pi/4, 3pi/4, -pi/4, -3pi/4.
    last = math.pi/4
    # Generate the total_phase 
    for i in range(len(phase_delta)):
        # Current phase is the current phase delta (phase change in this sample)
        # plus the absolute phase from the last sample
        a = phase_delta[i]+last
        # Normalize to keep between -pi and pi
        a = normalizeAngleRadians(a)
        total_phase.append(a)
        # Update the last phase value, need it for next cycle
        last = total_phase[i]

    # Write the modulated total_phase to a file.
    of = open('phase.txt', 'w')
    for i in total_phase:
        of.write(str(int(1024*math.cos(i))))
        of.write(' ')
        of.write(str(int(1024*math.sin(i))))
        of.write('\n')
    of.close()

##########################################################################
# SIMPLE RECEIVE PROCESSING - DEMODULATION
# Demodulation method used here is as follows:
# Divide every symbol interval into two bit intervals, intervals are 
# I bit and Q bit. The shaping filter will cause rotations at a bit 
# interval rate.
# The even bit interval is mapped onto the real axis and the odd bit interval
# is mapped onto the imaginary axis.
# When making an I bit decision, all that matters is the real value of the IQ 
# point. If greater than 0 (angle -pi/2 to pi/2) then call I for 0 else
# if less than 0 call I for 1
# When making a Q bit decision, all that matters is the imag value of the IQ 
# point. If greater than 0 (angle 0 to pi) then call Q for 0 else
# if less than 0 call Q for 1
##########################################################################

#                  IQ Constellation
#                        Q
#                        |
#                        |
#            10          |         00
#                        |
#                        |
#                        |
#                        |
#------------------------|------------------------ I
#                        |
#                        |
#                        |
#            11          |         01
#                        |
#                        |
#                        |
#                        |

#               Even Bit Processing - I Bit
#                        Q
#                        |
#                        |
#            1           |         0
#                        |
#                        |
#                        |
#                        |
#------------------------|------------------------ I
#                        |
#                        |
#                        |
#            1           |         0
#                        |
#                        |
#                        |
#                        |

#               Odd Bit Processing - Q Bit
#                        Q
#                        |
#                        |
#            0           |          0
#                        |
#                        |
#                        |
#                        |
#------------------------|------------------------ I
#                        |
#                        |
#                        |
#             1          |          1
#                        |
#                        |
#                        |
#                        |


    # Keep track of even/odd bit interval
    isEven = 1
    demod_bits = []
    # Demodulation is occurring here as described above
    # Incrementing along at the bit interval and making soft decisions based on
    # Quadrant that the IQ sample is located in and the bit interval
    # for loop paramaters:
    #################################################################################################
    # I BIT # Q BIT # I BIT # Q BIT # I BIT # Q BIT # I BIT # Q BIT # I BIT # Q BIT # I BIT # Q BIT #                                                                   
    #################################################################################################
    #   SYMBOL INT  #   SYMBOL INT  #   SYMBOL INT  #   SYMBOL INT  #   SYMBOL INT  #   SYMBOL INT  #
    #################################################################################################
    #   start = SAMPLES_PER_SYMBOL//4 = Want to capture each bit in the middle of bit interval.
    #                                   SAMPLES_PER_SYMBOL//2 would be right at transition of
    #                                   bit interval, so to capture middle of bit interval
    #                                   this needs to be SAMPLES_PER_SYMBOL//4
    #   end   = len(total_phase)      = Loop over all samples.. i.e. keep incrementing until all
    #                                   samples are processed
    #   step  = SAMPLES_PER_SYMBOL//2 = Want to step through samples bit by bit, So need to increment
    #                                   along at half the symbol rate to be at the bit rate
    for i in range(SAMPLES_PER_SYMBOL//4,len(total_phase), SAMPLES_PER_SYMBOL//2):
        # This is the total phase, aka the constellation point (angle) for the current IQ sample
        # to be used in the bit decision process
        a = total_phase[i]
        # If we are processing an I Bit
        if isEven == 1:
            # Bit decision threshold is the imag axis.
            # If we are to the right of the imag axis, call for 0
            if (a > -math.pi/2) and (a < math.pi/2):
                b = 0
            # If we are to the left of the imag axis, call for 1
            else:
                b = 1
        # If we are processing a Q bit
        else:
            # Bit decision threshold is the real axis.
            # If we are above the real axis, call for 0
            if (a < math.pi) and (a > 0):
                b = 0
            # If we are below the real axis, call for 1
            else:
                b = 1
        # Append the bits that have been demodulated to array of demodulated bit
        demod_bits.append(b)
        # Toggle the bit interval between even and odd
        isEven ^= 1

    # Need to take the first and last 6 bits of the end of the demodulated array
    # This is due to the nature of the shaping filter. The shaping filter convolution runs over 12 bits,
    # (although some are zero valued) which is 6 symbol periods. Therefore this convolution introduces
    # a 6 extra bits on both the end and start of the demodulated bits.
    demod_bits = demod_bits[6:-6]

    print("Transmitted Bits")
    #print(bits)
    print("Demodulated Bits")
    #print(demod_bits)
    if bits==demod_bits:
        print("TX bits match demodulated bits")
    else:
        print("TX bits DO NOT match demodulated bits")


if __name__ == "__main__":
    main()
