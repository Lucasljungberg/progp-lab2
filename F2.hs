--
-- Labb F2: Molykylärbiologi i Haskell
-- Progp 15 
-- Peter Jonsson <peterjo6@kth.se> och Lucas Ljungberg <lucaslj@kth.se>
-- KTH D-14
-- 2015-09-07
--

-- Required module F2
module F2 where
    import Data.Bits
    import Data.List
    import Data.Maybe

    data MolSeq = MolSeq {name :: String, mSeq :: String, seqType :: Type}
    data Type = DNA | PROTEIN deriving (Eq, Enum)

    string2seq :: String -> String -> MolSeq
    string2seq name seq = MolSeq name seq (checkType seq) 

    seqName :: MolSeq -> String
    seqName molSeq = (name molSeq)

    seqSequence :: MolSeq -> String
    seqSequence molSeq = (mSeq molSeq)

    seqLength :: MolSeq -> Int
    seqLength molSeq = length (mSeq molSeq)

    -- Calculates the evolutionary distance of two MolSeq instances.
    -- This is done using the Jukes-Cantor model for DNA sequences and the Poisson model for proteins.
    -- This process is O(n)
    seqDistance :: MolSeq -> MolSeq -> Double
    seqDistance mol1 mol2 = if (seqType mol1) /= (seqType mol2)
        then error("Cannot compare DNA with proteins") else d where
            d = if isProtein mol1 && a >= 0.94 then 3.7
                else
                    if a > 0.74 then 3.3
                    else -((r-1)/r) * log (1 - ((r * a)/(r-1)))
            r = if isProtein mol1 then 20 else 4
            a = fromIntegral(sum differences) / fromIntegral(length seq1)
            seq1 = (mSeq mol1)
            seq2 = (mSeq mol2)
            differences = [if seq1!!n /= seq2!!n then 1 else 0 | n <- [0..(length seq1)-1]]


    checkType :: String -> Type
    checkType seq = is
        where
            is = if length (filter (`notElem` "ACGT") seq) /= 0 then PROTEIN else DNA
                 
    isProtein :: MolSeq -> Bool
    isProtein seq = if (seqType seq) == PROTEIN then True else False

    data Profile = Profile {pName :: String, matrix :: [[(Char, Int)]], numSeq :: Int, profType :: Type}

    molseqs2profile :: String -> [MolSeq] -> Profile
    molseqs2profile name molSeq = Profile name matrix numSeq isP
        where
            numSeq = length molSeq
            matrix = makeProfileMatrix molSeq
            isP = checkType (mSeq (head molSeq))

    profileName :: Profile -> String
    profileName profile = pName profile

    profileFrequency :: Profile -> Int -> Char -> Double
    profileFrequency profile pos token = count / total
        where
            tokenPos = fromMaybe (error "Token not found") (elemIndex token list)
            list = if (profType profile) == PROTEIN then aminoacids else nucleotides
            innerList = ((matrix profile)!!pos)
            count = fromIntegral (snd (innerList!!tokenPos))
            total = fromIntegral (numSeq profile)

    -- Gets the evolutionary distance between two profiles.
    profileDistance :: Profile -> Profile -> Double
    profileDistance p1 p2 = sum distances where
        distances = [abs ((profileFrequency p1 (fromIntegral j) i) - (profileFrequency p2 (fromIntegral j) i)) | i <- list, j <- [0..stop]]
        list = if (profType p1) == PROTEIN then aminoacids else nucleotides
        stop = length (matrix p1) - 1


    nucleotides = "ACGT"
    aminoacids = sort "ARNDCEQGHILKMFPSTWYVX"

    makeProfileMatrix :: [MolSeq] -> [[(Char, Int)]]
    makeProfileMatrix [] = error "Empty sequence list"
    makeProfileMatrix sl = res
      where 
        t = seqType (head sl)
        defaults = 
          if (t == DNA) then
            zip nucleotides (replicate (length nucleotides) 0) -- Rad (i)
          else 
            zip aminoacids (replicate (length aminoacids) 0)   -- Rad (ii)
        strs = map seqSequence sl                              -- Rad (iii)
        tmp1 = map (map (\x -> ((head x), (length x))) . group . sort)
                   (transpose strs)                            -- Rad (iv)
        equalFst a b = (fst a) == (fst b)
        res = map sort (map (\l -> unionBy equalFst l defaults) tmp1) 
