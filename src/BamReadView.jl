
import Base: copy, start, next, done, eltype, print

# ViewTypes are just used a decorators
abstract ViewType
abstract ViewTypeUInt8 <: ViewType   # start in bytes, length in bytes
abstract ViewTypeUInt4 <: ViewType   # start in halfbytes, length in halfbytes
abstract ViewTypeUInt32 <: ViewType  # start in bytes, length in UInt32s


immutable SeqElement2
	e::Int
end
const base_indices = [0,1,2,0,3,0,0,0,4,0,0,0,0,0,0,0]
const base_chars = UInt8['=','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'];
base_mask(e::SeqElement2) = e.e
base_index(e::SeqElement2) = base_indices[e.e+1]
base_char(e::SeqElement2) = base_chars[e.e+1] # TODO: should this be renamed, return Char or change in some other way?
print(io::IO,e::SeqElement2) = print(io,Char(base_char(e)))


immutable CIGARElement2
	e::UInt32
end
op(e::CIGARElement2) = Int(e.e&0xf)
const op_chars = Char['M','I','D','N','S','H','P','=','X']
op_char(e::CIGARElement2) = op_chars[op(e)+1]
len(e::CIGARElement2) = Int(e.e>>4)
print(io::IO,e::CIGARElement2) = print(io,len(e),op_char(e))

# "concrete" view types
abstract ViewTypeSeq <: ViewTypeUInt4
abstract ViewTypeCIGAR <: ViewTypeUInt32



immutable BamReadView{T<:ViewType}
	buf::Array{UInt8,1} # entire buffer - typically the same as the BamRead buffer
	start::Int          # start offset  - interpretation depends on T (1-based index)
	len::Int            # length        - interpretation depends on T
end


typealias Seq2 BamReadView{ViewTypeSeq}
seq2(r::BamRead) = Seq2(r.buf, (r.seq_start-1)*2+1, l_seq(r))


typealias CIGAR2 BamReadView{ViewTypeCIGAR}
cigar2(r::BamRead) = CIGAR2(r.buf, l_read_name(r)+1, n_cigar_op(r))



copy{T<:ViewTypeUInt8}(v::BamReadView{T}) = BamReadView{T}( v.buf[v.start:v.start+v.len-1], 1, v.len )
function copy{T<:ViewTypeUInt4}(v::BamReadView{T})
	byterange = (v.start-1)>>1+1:(v.start+v.len)>>1
	BamReadView{T}( v.buf[byterange], mod(v.start-1,2)+1, v.len )
end
copy{T<:ViewTypeUInt32}(v::BamReadView{T}) = BamReadView{T}( v.buf[v.start:v.start+v.len*4-1], 1, v.len )



# iterator with helper function for interpretation of value

getvalue{T<:ViewTypeUInt8}(v::BamReadView{T},i::Int) = v.buf[i+v.start]
function getvalue{T<:ViewTypeUInt4}(v::BamReadView{T},i::Int)
	i2 = v.start+i+1
	b = v.buf[i2>>1]
	i2&1==0 ? b>>4 : b&0x0f # select low or high nybble depending on index
end
function getvalue{T<:ViewTypeUInt32}(v::BamReadView{T},i::Int)
	i2 = v.start+i*4
	UInt32(v.buf[i2]) | (UInt32(v.buf[i2+1])<<8) | (UInt32(v.buf[i2+2])<<16) | (UInt32(v.buf[i2+3])<<24)
end



interpret{T<:ViewType}(v::BamReadView{T},value) = value
interpret{T<:ViewTypeSeq}(v::BamReadView{T},value) = SeqElement2(value)
interpret{T<:ViewTypeCIGAR}(v::BamReadView{T},value) = CIGARElement2(value)


start{T<:ViewType}(::BamReadView{T}) = 0 # 0-based sinc v.start is 1-based
next{T<:ViewType}(v::BamReadView{T}, i::Int) = (interpret(v,getvalue(v,i)), i+1)
done{T<:ViewType}(s::BamReadView{T}, i::Int) = i >= s.len
#eltype(::Type{BamReadView{T}}) = 


eltype{T<:ViewTypeUInt8}(::Type{BamReadView{T}}) = UInt8
eltype{T<:ViewTypeUInt32}(::Type{BamReadView{T}}) = UInt32
eltype{T<:ViewTypeSeq}(::Type{BamReadView{T}}) = SeqElement2
eltype{T<:ViewTypeCIGAR}(::Type{BamReadView{T}}) = CIGARElement2




# array functionality
